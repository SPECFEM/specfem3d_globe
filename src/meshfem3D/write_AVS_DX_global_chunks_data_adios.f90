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

! create AVS or DX 2D data for the faces of the global chunks,
! to be recombined in postprocessing

module AVS_DX_global_chunks_mod

  implicit none

  type avs_dx_global_chunks_t
    integer(kind=4) :: npoin, nspecface
    real(kind=4), dimension(:), allocatable :: x_adios, y_adios, z_adios
    integer(kind=4), dimension(:), allocatable :: idoubling, iglob1, iglob2, &
        iglob3, iglob4
    real(kind=4), dimension(:), allocatable :: vmin, vmax
    real(kind=4), dimension(:), allocatable :: dvp, dvs
  endtype

contains


  subroutine define_AVS_DX_global_chunks_data(nspec,iboun,ibool, &
                                              mask_ibool, &
                                              npointot, &
                                              MODEL_3D_MANTLE_PERTUBATIONS, &
                                              group_size_inc, avs_dx_adios)

  use constants
  use adios_helpers_mod
  use manager_adios

  implicit none

  integer :: nspec

  logical :: iboun(6,nspec)
  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  integer :: npointot
  ! logical mask used to output global points only once
  logical :: mask_ibool(npointot)

  logical :: MODEL_3D_MANTLE_PERTUBATIONS

  integer(kind=8), intent(inout) :: group_size_inc
  type(avs_dx_global_chunks_t), intent(inout) :: avs_dx_adios

  ! local parameters
  integer :: ispec
  integer, dimension(8) :: iglobval
  integer :: npoin,nspecface
  integer(kind=8) :: local_dim

  integer :: ierr

  ! Dummy arrays for type inference inside adios helpers, need correct size as well
  real(kind=4), dimension(:), allocatable :: dummy_real1d
  integer(kind=4), dimension(:), allocatable :: dummy_int1d

  mask_ibool(:) = .false.

  nspecface = 0

  ! mark global AVS or DX points
  do ispec = 1,nspec
  ! only if on face
    if (iboun(1,ispec) .or. iboun(2,ispec) .or. &
                iboun(3,ispec) .or. iboun(4,ispec)) then
      iglobval(1)=ibool(1,1,1,ispec)
      iglobval(2)=ibool(NGLLX,1,1,ispec)
      iglobval(3)=ibool(NGLLX,NGLLY,1,ispec)
      iglobval(4)=ibool(1,NGLLY,1,ispec)
      iglobval(5)=ibool(1,1,NGLLZ,ispec)
      iglobval(6)=ibool(NGLLX,1,NGLLZ,ispec)
      iglobval(7)=ibool(NGLLX,NGLLY,NGLLZ,ispec)
      iglobval(8)=ibool(1,NGLLY,NGLLZ,ispec)

      ! face xi = xi_min
      if (iboun(1,ispec)) then
        nspecface = nspecface + 1
        mask_ibool(iglobval(1)) = .true.
        mask_ibool(iglobval(4)) = .true.
        mask_ibool(iglobval(8)) = .true.
        mask_ibool(iglobval(5)) = .true.
      endif

      ! face xi = xi_max
      if (iboun(2,ispec)) then
        nspecface = nspecface + 1
        mask_ibool(iglobval(2)) = .true.
        mask_ibool(iglobval(3)) = .true.
        mask_ibool(iglobval(7)) = .true.
        mask_ibool(iglobval(6)) = .true.
      endif

      ! face eta = eta_min
      if (iboun(3,ispec)) then
        nspecface = nspecface + 1
        mask_ibool(iglobval(1)) = .true.
        mask_ibool(iglobval(2)) = .true.
        mask_ibool(iglobval(6)) = .true.
        mask_ibool(iglobval(5)) = .true.
      endif

      ! face eta = eta_max
      if (iboun(4,ispec)) then
        nspecface = nspecface + 1
        mask_ibool(iglobval(4)) = .true.
        mask_ibool(iglobval(3)) = .true.
        mask_ibool(iglobval(7)) = .true.
        mask_ibool(iglobval(8)) = .true.
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

  allocate(avs_dx_adios%vmin(npoin), stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, "Error allocating vmin.")
  allocate(avs_dx_adios%vmax(npoin), stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, "Error allocating vmax.")

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

  !--- Variables for '...AVS_DXpointschunks.txt'
  local_dim = npoin
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, &
                                   '', "points_chunks/x_value", dummy_real1d)
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, &
                                   '', "points_chunks/y_value", dummy_real1d)
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, &
                                   '', "points_chunks/z_value", dummy_real1d)

  !--- Variables for '...AVS_DXpointschunk_stability.txt'
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, &
                                   '', "points_chunks/vmin", dummy_real1d)
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, &
                                   '', "points_chunks/vmax", dummy_real1d)

!  !--- Variables for AVS_DXelementschunks.txt
  local_dim = nspecface
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, &
                                   '', "elements_chunks/idoubling", dummy_int1d)

  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, &
                                   '', "elements_chunks/num_ibool_AVS_DX_iglob1", dummy_int1d)
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, &
                                   '', "elements_chunks/num_ibool_AVS_DX_iglob2", dummy_int1d)
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, &
                                   '', "elements_chunks/num_ibool_AVS_DX_iglob3", dummy_int1d)
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, &
                                   '', "elements_chunks/num_ibool_AVS_DX_iglob4", dummy_int1d)

  if (MODEL_3D_MANTLE_PERTUBATIONS) then
    allocate(avs_dx_adios%dvp(nspecface), stat=ierr)
    if (ierr /= 0) call exit_MPI(myrank, "Error allocating dvp.")
    allocate(avs_dx_adios%dvs(nspecface), stat=ierr)
    if (ierr /= 0) call exit_MPI(myrank, "Error allocating dvs.")

    ! re-allocates dummy for definitions to have correct size
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

  end subroutine define_AVS_DX_global_chunks_data

!===============================================================================

  subroutine prepare_AVS_DX_global_chunks_data_adios(prname,nspec, &
                                                     iboun,ibool, idoubling,xstore,ystore,zstore,num_ibool_AVS_DX,mask_ibool, &
                                                     npointot,rhostore,kappavstore,muvstore, &
                                                     nspl,rspl,ellipicity_spline,ellipicity_spline2,ELLIPTICITY, &
                                                     MODEL_3D_MANTLE_PERTUBATIONS, &
                                                     RICB,RCMB,RTOPDDOUBLEPRIME,R670,R220,R771,R400,R120,R80,RMOHO, &
                                                     RMIDDLE_CRUST,iregion_code, &
                                                     avs_dx_adios)

  use constants
  use shared_parameters, only: R_PLANET

  implicit none

  ! processor identification
  character(len=MAX_STRING_LEN) :: prname

  integer :: nspec

  logical iboun(6,nspec)

  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  integer idoubling(nspec)

  double precision,dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

  integer :: npointot
  ! numbering of global AVS or DX points
  integer num_ibool_AVS_DX(npointot)
  ! logical mask used to output global points only once
  logical mask_ibool(npointot)

  real(kind=CUSTOM_REAL) kappavstore(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) muvstore(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) rhostore(NGLLX,NGLLY,NGLLZ,nspec)

  ! for ellipticity
  integer nspl
  double precision rspl(NR_DENSITY),ellipicity_spline(NR_DENSITY),ellipicity_spline2(NR_DENSITY)

  logical ELLIPTICITY,MODEL_3D_MANTLE_PERTUBATIONS

  double precision RICB,RCMB,RTOPDDOUBLEPRIME,R670,R220,R771, &
    R400,R120,R80,RMOHO,RMIDDLE_CRUST

  integer iregion_code

  ! local parameters
  integer ispec
  integer i,j,k,np
  integer, dimension(8) :: iglobval
  integer npoin,numpoin,nspecface,ispecface

  real(kind=CUSTOM_REAL) vmin,vmax

  double precision r,rho,vp,vs,Qkappa,Qmu
  double precision vpv,vph,vsv,vsh,eta_aniso
  double precision x,y,z
  real(kind=CUSTOM_REAL) dvp,dvs

  type(avs_dx_global_chunks_t), intent(inout) :: avs_dx_adios ! out for adios_write

  integer :: idummy

  ! to avoid compiler warnings
  idummy = len_trim(prname)

  ! erase the logical mask used to mark points already found
  mask_ibool(:) = .false.

  nspecface = 0
  dvp = 0.0
  dvs = 0.0

  ! mark global AVS or DX points
  do ispec = 1,nspec
  ! only if on face
    if (iboun(1,ispec) .or. iboun(2,ispec) .or. &
                iboun(3,ispec) .or. iboun(4,ispec)) then
      iglobval(1)=ibool(1,1,1,ispec)
      iglobval(2)=ibool(NGLLX,1,1,ispec)
      iglobval(3)=ibool(NGLLX,NGLLY,1,ispec)
      iglobval(4)=ibool(1,NGLLY,1,ispec)
      iglobval(5)=ibool(1,1,NGLLZ,ispec)
      iglobval(6)=ibool(NGLLX,1,NGLLZ,ispec)
      iglobval(7)=ibool(NGLLX,NGLLY,NGLLZ,ispec)
      iglobval(8)=ibool(1,NGLLY,NGLLZ,ispec)

      ! face xi = xi_min
      if (iboun(1,ispec)) then
        nspecface = nspecface + 1
        mask_ibool(iglobval(1)) = .true.
        mask_ibool(iglobval(4)) = .true.
        mask_ibool(iglobval(8)) = .true.
        mask_ibool(iglobval(5)) = .true.
      endif

      ! face xi = xi_max
      if (iboun(2,ispec)) then
        nspecface = nspecface + 1
        mask_ibool(iglobval(2)) = .true.
        mask_ibool(iglobval(3)) = .true.
        mask_ibool(iglobval(7)) = .true.
        mask_ibool(iglobval(6)) = .true.
      endif

      ! face eta = eta_min
      if (iboun(3,ispec)) then
        nspecface = nspecface + 1
        mask_ibool(iglobval(1)) = .true.
        mask_ibool(iglobval(2)) = .true.
        mask_ibool(iglobval(6)) = .true.
        mask_ibool(iglobval(5)) = .true.
      endif

      ! face eta = eta_max
      if (iboun(4,ispec)) then
        nspecface = nspecface + 1
        mask_ibool(iglobval(4)) = .true.
        mask_ibool(iglobval(3)) = .true.
        mask_ibool(iglobval(7)) = .true.
        mask_ibool(iglobval(8)) = .true.
      endif

    endif
  enddo

  ! count global number of AVS or DX points
  npoin = count(mask_ibool(:))

  ! safety check
  if (npoin /= avs_dx_adios%npoin ) stop 'Error npoin invalid in prepare_AVS_DX_global_chunks_data_adios() routine'

  ! number of points in AVS or DX file
  !write(*,*) npoin

  ! erase the logical mask used to mark points already found
  mask_ibool(:) = .false.

  ! output global AVS or DX points
  numpoin = 0
  do ispec = 1,nspec
  ! only if on face
    if (iboun(1,ispec) .or. iboun(2,ispec) .or. &
                iboun(3,ispec) .or. iboun(4,ispec)) then
      iglobval(1)=ibool(1,1,1,ispec)
      iglobval(2)=ibool(NGLLX,1,1,ispec)
      iglobval(3)=ibool(NGLLX,NGLLY,1,ispec)
      iglobval(4)=ibool(1,NGLLY,1,ispec)
      iglobval(5)=ibool(1,1,NGLLZ,ispec)
      iglobval(6)=ibool(NGLLX,1,NGLLZ,ispec)
      iglobval(7)=ibool(NGLLX,NGLLY,NGLLZ,ispec)
      iglobval(8)=ibool(1,NGLLY,NGLLZ,ispec)

      ! face xi = xi_min
      if (iboun(1,ispec)) then

        if (.not. mask_ibool(iglobval(1))) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglobval(1)) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(1,1,1,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(1,1,1,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(1,1,1,ispec))

          vmax = sqrt((kappavstore(1,1,1,ispec) &
                       + 4.*muvstore(1,1,1,ispec)/3.)/rhostore(1,1,1,ispec))
          vmin = sqrt(muvstore(1,1,1,ispec)/rhostore(1,1,1,ispec))
          ! particular case of the outer core (muvstore contains 1/rho)
          if (idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL) then
            r = dsqrt(xstore(1,1,1,ispec)**2 + ystore(1,1,1,ispec)**2 &
                + zstore(1,1,1,ispec)**2)
            call prem_display_outer_core(r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
            vmax = real(vp,kind=CUSTOM_REAL)
            vmin = real(vp,kind=CUSTOM_REAL)
          endif
          if (vmin == 0.0) vmin = vmax
          avs_dx_adios%vmin(numpoin) = vmin
          avs_dx_adios%vmax(numpoin) = vmax
        endif

        if (.not. mask_ibool(iglobval(4))) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglobval(4)) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(1,NGLLY,1,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(1,NGLLY,1,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(1,NGLLY,1,ispec))

          vmax = sqrt((kappavstore(1,NGLLY,1,ispec) &
                       +4.*muvstore(1,NGLLY,1,ispec)/3.)/rhostore(1,NGLLY,1,ispec))
          vmin = sqrt(muvstore(1,NGLLY,1,ispec)/rhostore(1,NGLLY,1,ispec))
          ! particular case of the outer core (muvstore contains 1/rho)
          if (idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL) then
            r = dsqrt(xstore(1,NGLLY,1,ispec)**2 + ystore(1,NGLLY,1,ispec)**2 &
                + zstore(1,NGLLY,1,ispec)**2)
            call prem_display_outer_core(r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
            vmax = real(vp,kind=CUSTOM_REAL)
            vmin = real(vp,kind=CUSTOM_REAL)
          endif
          if (vmin == 0.0) vmin = vmax
          avs_dx_adios%vmin(numpoin) = vmin
          avs_dx_adios%vmax(numpoin) = vmax
        endif

        if (.not. mask_ibool(iglobval(8))) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglobval(8)) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(1,NGLLY,NGLLZ,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(1,NGLLY,NGLLZ,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(1,NGLLY,NGLLZ,ispec))

          vmax = sqrt((kappavstore(1,NGLLY,NGLLZ,ispec) &
                       +4.*muvstore(1,NGLLY,NGLLZ,ispec)/3.) / rhostore(1,NGLLY,NGLLZ,ispec))
          vmin = sqrt(muvstore(1,NGLLY,NGLLZ,ispec) / rhostore(1,NGLLY,NGLLZ,ispec))
          ! particular case of the outer core (muvstore contains 1/rho)
          if (idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL) then
            r = dsqrt(xstore(1,NGLLY,NGLLZ,ispec)**2 &
                + ystore(1,NGLLY,NGLLZ,ispec)**2 &
                + zstore(1,NGLLY,NGLLZ,ispec)**2)
            call prem_display_outer_core(r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
            vmax = real(vp,kind=CUSTOM_REAL)
            vmin = real(vp,kind=CUSTOM_REAL)
          endif
          if (vmin == 0.0) vmin = vmax

          avs_dx_adios%vmin(numpoin) = vmin
          avs_dx_adios%vmax(numpoin) = vmax
        endif

        if (.not. mask_ibool(iglobval(5))) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglobval(5)) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(1,1,NGLLZ,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(1,1,NGLLZ,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(1,1,NGLLZ,ispec))

          vmax = sqrt((kappavstore(1,1,NGLLZ,ispec) &
                       +4.*muvstore(1,1,NGLLZ,ispec)/3.)/rhostore(1,1,NGLLZ,ispec))
          vmin = sqrt(muvstore(1,1,NGLLZ,ispec)/rhostore(1,1,NGLLZ,ispec))
          ! particular case of the outer core (muvstore contains 1/rho)
          if (idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL) then
            r = dsqrt(xstore(1,1,NGLLZ,ispec)**2 + ystore(1,1,NGLLZ,ispec)**2 &
                + zstore(1,1,NGLLZ,ispec)**2)
            call prem_display_outer_core(r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
            vmax = real(vp,kind=CUSTOM_REAL)
            vmin = real(vp,kind=CUSTOM_REAL)
          endif
          if (vmin == 0.0) vmin = vmax
          avs_dx_adios%vmin(numpoin) = vmin
          avs_dx_adios%vmax(numpoin) = vmax
        endif

        mask_ibool(iglobval(1)) = .true.
        mask_ibool(iglobval(4)) = .true.
        mask_ibool(iglobval(8)) = .true.
        mask_ibool(iglobval(5)) = .true.
      endif

      ! face xi = xi_max
      if (iboun(2,ispec)) then

        if (.not. mask_ibool(iglobval(2))) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglobval(2)) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(NGLLX,1,1,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(NGLLX,1,1,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(NGLLX,1,1,ispec))

          vmax = sqrt((kappavstore(NGLLX,1,1,ispec) &
                       +4.*muvstore(NGLLX,1,1,ispec)/3.)/rhostore(NGLLX,1,1,ispec))
          vmin = sqrt(muvstore(NGLLX,1,1,ispec)/rhostore(NGLLX,1,1,ispec))
          ! particular case of the outer core (muvstore contains 1/rho)
          if (idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL) then
            r = dsqrt(xstore(NGLLX,1,1,ispec)**2 + ystore(NGLLX,1,1,ispec)**2 &
                + zstore(NGLLX,1,1,ispec)**2)
            call prem_display_outer_core(r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
            vmax = real(vp,kind=CUSTOM_REAL)
            vmin = real(vp,kind=CUSTOM_REAL)
          endif
          if (vmin == 0.0) vmin = vmax
          avs_dx_adios%vmin(numpoin) = vmin
          avs_dx_adios%vmax(numpoin) = vmax
        endif

        if (.not. mask_ibool(iglobval(3))) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglobval(3)) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(NGLLX,NGLLY,1,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(NGLLX,NGLLY,1,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(NGLLX,NGLLY,1,ispec))

          vmax = sqrt((kappavstore(NGLLX,NGLLY,1,ispec) &
                       + 4.*muvstore(NGLLX,NGLLY,1,ispec)/3.) / rhostore(NGLLX,NGLLY,1,ispec))
          vmin = sqrt(muvstore(NGLLX,NGLLY,1,ispec) / rhostore(NGLLX,NGLLY,1,ispec))
          ! particular case of the outer core (muvstore contains 1/rho)
          if (idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL) then
            r = dsqrt(xstore(NGLLX,NGLLY,1,ispec)**2 &
                + ystore(NGLLX,NGLLY,1,ispec)**2 &
                + zstore(NGLLX,NGLLY,1,ispec)**2)
            call prem_display_outer_core(r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
            vmax = real(vp,kind=CUSTOM_REAL)
            vmin = real(vp,kind=CUSTOM_REAL)
          endif
          if (vmin == 0.0) vmin = vmax
          avs_dx_adios%vmin(numpoin) = vmin
          avs_dx_adios%vmax(numpoin) = vmax
        endif

        if (.not. mask_ibool(iglobval(7))) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglobval(7)) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(NGLLX,NGLLY,NGLLZ,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(NGLLX,NGLLY,NGLLZ,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(NGLLX,NGLLY,NGLLZ,ispec))

          vmax = sqrt((kappavstore(NGLLX,NGLLY,NGLLZ,ispec) &
                       + 4.*muvstore(NGLLX,NGLLY,NGLLZ,ispec)/3.) / rhostore(NGLLX,NGLLY,NGLLZ,ispec))
          vmin = sqrt(muvstore(NGLLX,NGLLY,NGLLZ,ispec) / rhostore(NGLLX,NGLLY,NGLLZ,ispec))
          ! particular case of the outer core (muvstore contains 1/rho)
          if (idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL) then
            r = dsqrt(xstore(NGLLX,NGLLY,NGLLZ,ispec)**2 &
                + ystore(NGLLX,NGLLY,NGLLZ,ispec)**2 &
                + zstore(NGLLX,NGLLY,NGLLZ,ispec)**2)
            call prem_display_outer_core(r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
            vmax = real(vp,kind=CUSTOM_REAL)
            vmin = real(vp,kind=CUSTOM_REAL)
          endif
          if (vmin == 0.0) vmin = vmax
          avs_dx_adios%vmin(numpoin) = vmin
          avs_dx_adios%vmax(numpoin) = vmax
        endif

        if (.not. mask_ibool(iglobval(6))) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglobval(6)) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(NGLLX,1,NGLLZ,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(NGLLX,1,NGLLZ,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(NGLLX,1,NGLLZ,ispec))

          vmax = sqrt((kappavstore(NGLLX,1,NGLLZ,ispec) &
                       + 4.*muvstore(NGLLX,1,NGLLZ,ispec)/3.) / rhostore(NGLLX,1,NGLLZ,ispec))
          vmin = sqrt(muvstore(NGLLX,1,NGLLZ,ispec) / rhostore(NGLLX,1,NGLLZ,ispec))
          ! particular case of the outer core (muvstore contains 1/rho)
          if (idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL) then
            r = dsqrt(xstore(NGLLX,1,NGLLZ,ispec)**2 &
                + ystore(NGLLX,1,NGLLZ,ispec)**2 &
                + zstore(NGLLX,1,NGLLZ,ispec)**2)
            call prem_display_outer_core(r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
            vmax = real(vp,kind=CUSTOM_REAL)
            vmin = real(vp,kind=CUSTOM_REAL)
          endif
          if (vmin == 0.0) vmin = vmax
          avs_dx_adios%vmin(numpoin) = vmin
          avs_dx_adios%vmax(numpoin) = vmax
        endif

        mask_ibool(iglobval(2)) = .true.
        mask_ibool(iglobval(3)) = .true.
        mask_ibool(iglobval(7)) = .true.
        mask_ibool(iglobval(6)) = .true.
      endif

      ! face eta = eta_min
      if (iboun(3,ispec)) then

        if (.not. mask_ibool(iglobval(1))) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglobval(1)) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(1,1,1,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(1,1,1,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(1,1,1,ispec))

          vmax = sqrt((kappavstore(1,1,1,ispec) &
                       + 4.*muvstore(1,1,1,ispec)/3.)/rhostore(1,1,1,ispec))
          vmin = sqrt(muvstore(1,1,1,ispec)/rhostore(1,1,1,ispec))
          ! particular case of the outer core (muvstore contains 1/rho)
          if (idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL) then
            r = dsqrt(xstore(1,1,1,ispec)**2 &
                + ystore(1,1,1,ispec)**2 + zstore(1,1,1,ispec)**2)
            call prem_display_outer_core(r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
            vmax = real(vp,kind=CUSTOM_REAL)
            vmin = real(vp,kind=CUSTOM_REAL)
          endif
          if (vmin == 0.0) vmin = vmax
          avs_dx_adios%vmin(numpoin) = vmin
          avs_dx_adios%vmax(numpoin) = vmax
        endif

        if (.not. mask_ibool(iglobval(2))) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglobval(2)) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(NGLLX,1,1,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(NGLLX,1,1,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(NGLLX,1,1,ispec))

          vmax = sqrt((kappavstore(NGLLX,1,1,ispec) &
                       +4.*muvstore(NGLLX,1,1,ispec)/3.)/rhostore(NGLLX,1,1,ispec))
          vmin = sqrt(muvstore(NGLLX,1,1,ispec)/rhostore(NGLLX,1,1,ispec))
          ! particular case of the outer core (muvstore contains 1/rho)
          if (idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL) then
            r = dsqrt(xstore(NGLLX,1,1,ispec)**2 &
                + ystore(NGLLX,1,1,ispec)**2 + zstore(NGLLX,1,1,ispec)**2)
            call prem_display_outer_core(r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
            vmax = real(vp,kind=CUSTOM_REAL)
            vmin = real(vp,kind=CUSTOM_REAL)
          endif
          if (vmin == 0.0) vmin = vmax
          avs_dx_adios%vmin = vmin
          avs_dx_adios%vmax = vmax
        endif

        if (.not. mask_ibool(iglobval(6))) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglobval(6)) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(NGLLX,1,NGLLZ,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(NGLLX,1,NGLLZ,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(NGLLX,1,NGLLZ,ispec))

          vmax = sqrt((kappavstore(NGLLX,1,NGLLZ,ispec) &
                       + 4.*muvstore(NGLLX,1,NGLLZ,ispec)/3.) / rhostore(NGLLX,1,NGLLZ,ispec))
          vmin = sqrt(muvstore(NGLLX,1,NGLLZ,ispec) / rhostore(NGLLX,1,NGLLZ,ispec))
          ! particular case of the outer core (muvstore contains 1/rho)
          if (idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL) then
            r = dsqrt(xstore(NGLLX,1,NGLLZ,ispec)**2 &
                + ystore(NGLLX,1,NGLLZ,ispec)**2 &
                + zstore(NGLLX,1,NGLLZ,ispec)**2)
            call prem_display_outer_core(r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
            vmax = real(vp,kind=CUSTOM_REAL)
            vmin = real(vp,kind=CUSTOM_REAL)
          endif
          if (vmin == 0.0) vmin = vmax
          avs_dx_adios%vmin(numpoin) = vmin
          avs_dx_adios%vmax(numpoin) = vmax
        endif

        if (.not. mask_ibool(iglobval(5))) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglobval(5)) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(1,1,NGLLZ,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(1,1,NGLLZ,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(1,1,NGLLZ,ispec))

          vmax = sqrt((kappavstore(1,1,NGLLZ,ispec) &
                       + 4.*muvstore(1,1,NGLLZ,ispec)/3.) / rhostore(1,1,NGLLZ,ispec))
          vmin = sqrt(muvstore(1,1,NGLLZ,ispec)/rhostore(1,1,NGLLZ,ispec))
          ! particular case of the outer core (muvstore contains 1/rho)
          if (idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL) then
            r = dsqrt(xstore(1,1,NGLLZ,ispec)**2 &
                + ystore(1,1,NGLLZ,ispec)**2 + zstore(1,1,NGLLZ,ispec)**2)
            call prem_display_outer_core(r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
            vmax = real(vp,kind=CUSTOM_REAL)
            vmin = real(vp,kind=CUSTOM_REAL)
          endif
          if (vmin == 0.0) vmin = vmax
          avs_dx_adios%vmin(numpoin) = vmin
          avs_dx_adios%vmax(numpoin) = vmax
        endif

        mask_ibool(iglobval(1)) = .true.
        mask_ibool(iglobval(2)) = .true.
        mask_ibool(iglobval(6)) = .true.
        mask_ibool(iglobval(5)) = .true.
      endif

      ! face eta = eta_max
      if (iboun(4,ispec)) then

        if (.not. mask_ibool(iglobval(4))) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglobval(4)) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(1,NGLLY,1,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(1,NGLLY,1,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(1,NGLLY,1,ispec))

          vmax = sqrt((kappavstore(1,NGLLY,1,ispec) &
                       + 4.*muvstore(1,NGLLY,1,ispec)/3.)/rhostore(1,NGLLY,1,ispec))
          vmin = sqrt(muvstore(1,NGLLY,1,ispec)/rhostore(1,NGLLY,1,ispec))
          ! particular case of the outer core (muvstore contains 1/rho)
          if (idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL) then
            r = dsqrt(xstore(1,NGLLY,1,ispec)**2 &
                + ystore(1,NGLLY,1,ispec)**2 + zstore(1,NGLLY,1,ispec)**2)
            call prem_display_outer_core(r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
            vmax = real(vp,kind=CUSTOM_REAL)
            vmin = real(vp,kind=CUSTOM_REAL)
          endif
          if (vmin == 0.0) vmin = vmax
          avs_dx_adios%vmin(numpoin) = vmin
          avs_dx_adios%vmax(numpoin) = vmax
        endif

        if (.not. mask_ibool(iglobval(3))) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglobval(3)) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(NGLLX,NGLLY,1,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(NGLLX,NGLLY,1,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(NGLLX,NGLLY,1,ispec))

          vmax = sqrt((kappavstore(NGLLX,NGLLY,1,ispec) &
                       + 4.*muvstore(NGLLX,NGLLY,1,ispec)/3.) / rhostore(NGLLX,NGLLY,1,ispec))
          vmin = sqrt(muvstore(NGLLX,NGLLY,1,ispec) / rhostore(NGLLX,NGLLY,1,ispec))
          ! particular case of the outer core (muvstore contains 1/rho)
          if (idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL) then
            r = dsqrt(xstore(NGLLX,NGLLY,1,ispec)**2 &
                + ystore(NGLLX,NGLLY,1,ispec)**2 &
                + zstore(NGLLX,NGLLY,1,ispec)**2)
            call prem_display_outer_core(r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
            vmax = real(vp,kind=CUSTOM_REAL)
            vmin = real(vp,kind=CUSTOM_REAL)
          endif

          if (vmin == 0.0) vmin = vmax

          avs_dx_adios%vmin(numpoin) = vmin
          avs_dx_adios%vmax(numpoin) = vmax
        endif

        if (.not. mask_ibool(iglobval(7))) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglobval(7)) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(NGLLX,NGLLY,NGLLZ,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(NGLLX,NGLLY,NGLLZ,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(NGLLX,NGLLY,NGLLZ,ispec))

          vmax = sqrt((kappavstore(NGLLX,NGLLY,NGLLZ,ispec) &
                       + 4.*muvstore(NGLLX,NGLLY,NGLLZ,ispec)/3.) / rhostore(NGLLX,NGLLY,NGLLZ,ispec))
          vmin = sqrt(muvstore(NGLLX,NGLLY,NGLLZ,ispec) / rhostore(NGLLX,NGLLY,NGLLZ,ispec))
          ! particular case of the outer core (muvstore contains 1/rho)
          if (idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL) then
            r = dsqrt(xstore(NGLLX,NGLLY,NGLLZ,ispec)**2 &
                + ystore(NGLLX,NGLLY,NGLLZ,ispec)**2 &
                + zstore(NGLLX,NGLLY,NGLLZ,ispec)**2)
            call prem_display_outer_core(r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
            vmax = real(vp,kind=CUSTOM_REAL)
            vmin = real(vp,kind=CUSTOM_REAL)
          endif
          if (vmin == 0.0) vmin = vmax
          avs_dx_adios%vmin(numpoin) = vmin
          avs_dx_adios%vmax(numpoin) = vmax
        endif

        if (.not. mask_ibool(iglobval(8))) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglobval(8)) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(1,NGLLY,NGLLZ,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(1,NGLLY,NGLLZ,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(1,NGLLY,NGLLZ,ispec))

          vmax = sqrt((kappavstore(1,NGLLY,NGLLZ,ispec) &
                       + 4.*muvstore(1,NGLLY,NGLLZ,ispec)/3.) / rhostore(1,NGLLY,NGLLZ,ispec))
          vmin = sqrt(muvstore(1,NGLLY,NGLLZ,ispec) / rhostore(1,NGLLY,NGLLZ,ispec))
          ! particular case of the outer core (muvstore contains 1/rho)
          if (idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL) then
            r = dsqrt(xstore(1,NGLLY,NGLLZ,ispec)**2 &
                + ystore(1,NGLLY,NGLLZ,ispec)**2 &
                + zstore(1,NGLLY,NGLLZ,ispec)**2)
            call prem_display_outer_core(r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
            vmax = real(vp,kind=CUSTOM_REAL)
            vmin = real(vp,kind=CUSTOM_REAL)
          endif
          if (vmin == 0.0) vmin = vmax

          avs_dx_adios%vmin(numpoin) = vmin
          avs_dx_adios%vmax(numpoin) = vmax
        endif

        mask_ibool(iglobval(4)) = .true.
        mask_ibool(iglobval(3)) = .true.
        mask_ibool(iglobval(7)) = .true.
        mask_ibool(iglobval(8)) = .true.
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
    ! only if on face
    if (iboun(1,ispec) .or. iboun(2,ispec) .or. &
        iboun(3,ispec) .or. iboun(4,ispec)) then
      iglobval(1)=ibool(1,1,1,ispec)
      iglobval(2)=ibool(NGLLX,1,1,ispec)
      iglobval(3)=ibool(NGLLX,NGLLY,1,ispec)
      iglobval(4)=ibool(1,NGLLY,1,ispec)
      iglobval(5)=ibool(1,1,NGLLZ,ispec)
      iglobval(6)=ibool(NGLLX,1,NGLLZ,ispec)
      iglobval(7)=ibool(NGLLX,NGLLY,NGLLZ,ispec)
      iglobval(8)=ibool(1,NGLLY,NGLLZ,ispec)

      ! include lateral variations if needed

      if (MODEL_3D_MANTLE_PERTUBATIONS) then
        !   pick a point within the element and get its radius
        r = dsqrt(xstore(2,2,2,ispec)**2+ystore(2,2,2,ispec)**2 &
            +zstore(2,2,2,ispec)**2)

        if (r > RCMB/R_PLANET .and. r < R_UNIT_SPHERE) then
          !     average over the element
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

                ! get reference model values: rho,vpv,vph,vsv,vsh and eta_aniso
                call meshfem3D_models_get1D_val(iregion_code, &
                                                idoubling(ispec), &
                                                r,rho,vpv,vph,vsv,vsh,eta_aniso, &
                                                Qkappa,Qmu,RICB,RCMB, &
                                                RTOPDDOUBLEPRIME,R80,R120,R220,R400,R670,R771, &
                                                RMOHO,RMIDDLE_CRUST)

                ! calculates isotropic values
                vp = sqrt(((8.d0+4.d0*eta_aniso)*vph*vph + 3.d0*vpv*vpv &
                        + (8.d0 - 8.d0*eta_aniso)*vsv*vsv)/15.d0)
                vs = sqrt(((1.d0-2.d0*eta_aniso)*vph*vph + vpv*vpv &
                        + 5.d0*vsh*vsh + (6.d0+4.d0*eta_aniso)*vsv*vsv)/15.d0)

                if (abs(rhostore(i,j,k,ispec)) < 1.e-20) then
                  print *,' attention: rhostore close to zero', &
                      rhostore(i,j,k,ispec),r,i,j,k,ispec
                  dvp = 0.0
                  dvs = 0.0
                else if (abs(sngl(vp)) < 1.e-20) then
                  print *,' attention: vp close to zero', &
                      sngl(vp),r,i,j,k,ispec
                  dvp = 0.0
                else if (abs(sngl(vs)) < 1.e-20) then
                  print *,' attention: vs close to zero', &
                      sngl(vs),r,i,j,k,ispec
                  dvs = 0.0
                else
                  dvp = dvp + (sqrt((kappavstore(i,j,k,ispec) &
                      +4.*muvstore(i,j,k,ispec)/3.)/rhostore(i,j,k,ispec)) &
                      - sngl(vp))/sngl(vp)
                  dvs = dvs + (sqrt(muvstore(i,j,k,ispec) &
                               / rhostore(i,j,k,ispec)) &
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
      if (iboun(1,ispec)) then
        ispecface = ispecface + 1
        avs_dx_adios%idoubling(ispecface) = idoubling(ispec)
        avs_dx_adios%iglob1(ispecface) = num_ibool_AVS_DX(iglobval(1))
        avs_dx_adios%iglob2(ispecface) = num_ibool_AVS_DX(iglobval(4))
        avs_dx_adios%iglob3(ispecface) = num_ibool_AVS_DX(iglobval(8))
        avs_dx_adios%iglob4(ispecface) = num_ibool_AVS_DX(iglobval(5))
        if (MODEL_3D_MANTLE_PERTUBATIONS) then
          avs_dx_adios%dvp(ispecface) = dvp
          avs_dx_adios%dvs(ispecface) = dvs
        endif
      endif

      ! face xi = xi_max
      if (iboun(2,ispec)) then
        ispecface = ispecface + 1
        avs_dx_adios%idoubling(ispecface) = idoubling(ispec)
        avs_dx_adios%iglob1(ispecface)= num_ibool_AVS_DX(iglobval(2))
        avs_dx_adios%iglob2(ispecface) = num_ibool_AVS_DX(iglobval(3))
        avs_dx_adios%iglob3(ispecface) = num_ibool_AVS_DX(iglobval(7))
        avs_dx_adios%iglob4(ispecface) = num_ibool_AVS_DX(iglobval(6))
        if (MODEL_3D_MANTLE_PERTUBATIONS) then
          avs_dx_adios%dvp(ispecface) = dvp
          avs_dx_adios%dvs(ispecface) = dvs
        endif
      endif

    ! face eta = eta_min
      if (iboun(3,ispec)) then
        ispecface = ispecface + 1
        avs_dx_adios%idoubling(ispecface) = idoubling(ispec)
        avs_dx_adios%iglob1(ispecface) = num_ibool_AVS_DX(iglobval(1))
        avs_dx_adios%iglob2(ispecface) = num_ibool_AVS_DX(iglobval(2))
        avs_dx_adios%iglob3(ispecface) = num_ibool_AVS_DX(iglobval(6))
        avs_dx_adios%iglob4(ispecface) = num_ibool_AVS_DX(iglobval(5))
        if (MODEL_3D_MANTLE_PERTUBATIONS) then
          avs_dx_adios%dvp(ispecface) = dvp
          avs_dx_adios%dvs(ispecface) = dvs
        endif
      endif

      ! face eta = eta_max
      if (iboun(4,ispec)) then
        ispecface = ispecface + 1
        avs_dx_adios%idoubling(ispecface) = idoubling(ispec)
        avs_dx_adios%iglob1(ispecface) = num_ibool_AVS_DX(iglobval(4))
        avs_dx_adios%iglob2(ispecface) = num_ibool_AVS_DX(iglobval(3))
        avs_dx_adios%iglob3(ispecface) = num_ibool_AVS_DX(iglobval(7))
        avs_dx_adios%iglob4(ispecface) = num_ibool_AVS_DX(iglobval(8))
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

  end subroutine prepare_AVS_DX_global_chunks_data_adios

!===============================================================================

  subroutine write_AVS_DX_global_chunks_data_adios(myrank, sizeprocs, avs_dx_adios, MODEL_3D_MANTLE_PERTUBATIONS)

  use adios_helpers_mod
  use manager_adios

  implicit none
  !--- Arguments
  integer, intent(in) :: myrank, sizeprocs
  type(avs_dx_global_chunks_t), intent(inout) :: avs_dx_adios ! out for adios_write
  logical MODEL_3D_MANTLE_PERTUBATIONS
  !--- Variables
  integer :: npoin, nspec
  integer(kind=8) :: local_dim

  npoin = avs_dx_adios%npoin
  nspec = avs_dx_adios%nspecface

  local_dim = npoin
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, &
                                   "points_chunks/x_value", avs_dx_adios%x_adios)
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, &
                                   "points_chunks/y_value", avs_dx_adios%y_adios)
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, &
                                   "points_chunks/z_value", avs_dx_adios%z_adios)

  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, &
                                   "points_chunks/vmin", avs_dx_adios%vmin)
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, &
                                   "points_chunks/vmax", avs_dx_adios%vmax)

  local_dim = nspec
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, &
                                  "elements_chunks/idoubling", avs_dx_adios%idoubling)

  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, &
                                   "elements_chunks/num_ibool_AVS_DX_iglob1", avs_dx_adios%iglob1)
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, &
                                   "elements_chunks/num_ibool_AVS_DX_iglob2", avs_dx_adios%iglob2)
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, &
                                   "elements_chunks/num_ibool_AVS_DX_iglob3", avs_dx_adios%iglob3)
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, &
                                   "elements_chunks/num_ibool_AVS_DX_iglob4", avs_dx_adios%iglob4)

  if (MODEL_3D_MANTLE_PERTUBATIONS) then
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, &
                                     "elements_faces/dvp", avs_dx_adios%dvp)
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, &
                                     "elements_faces/dvs", avs_dx_adios%dvs)
  endif

  end subroutine write_AVS_DX_global_chunks_data_adios

!===============================================================================

  subroutine free_AVS_DX_global_chunks_data_adios(avs_dx_adios, MODEL_3D_MANTLE_PERTUBATIONS)
  implicit none
  !--- Arguments
  type(avs_dx_global_chunks_t), intent(inout) :: avs_dx_adios
  logical MODEL_3D_MANTLE_PERTUBATIONS

  deallocate(avs_dx_adios%x_adios)
  deallocate(avs_dx_adios%y_adios)
  deallocate(avs_dx_adios%z_adios)

  deallocate(avs_dx_adios%vmin)
  deallocate(avs_dx_adios%vmax)

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

  end subroutine free_AVS_DX_global_chunks_data_adios

end module
