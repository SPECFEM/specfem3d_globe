!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
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

  subroutine save_arrays_solver(myrank,nspec,nglob,idoubling,ibool, &
                    iregion_code,xstore,ystore,zstore, &
                    NSPEC2D_TOP,NSPEC2D_BOTTOM)

  use constants

  use meshfem3D_models_par,only: &
    OCEANS,TRANSVERSE_ISOTROPY,HETEROGEN_3D_MANTLE,ANISOTROPIC_3D_MANTLE, &
    ANISOTROPIC_INNER_CORE,ATTENUATION

  use meshfem3D_par,only: &
    NCHUNKS,ABSORBING_CONDITIONS,SAVE_MESH_FILES

  use create_regions_mesh_par2,only: &
    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
    gammaxstore,gammaystore,gammazstore, &
    rhostore,dvpstore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
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
    prname

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
  integer :: NSPEC2D_TOP,NSPEC2D_BOTTOM

  ! local parameters
  integer :: i,j,k,ispec,iglob,ier
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: tmp_array

  ! mesh databases
  open(unit=27,file=prname(1:len_trim(prname))//'solver_data.bin', &
       status='unknown',form='unformatted',action='write',iostat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error opening solver_data.bin file')

  ! save nspec and nglob, to be used in combine_paraview_data
  write(27) nspec
  write(27) nglob

  ! mesh topology

  ! mesh arrays used in the solver to locate source and receivers
  ! and for anisotropy and gravity, save in single precision
  ! use tmp_array for temporary storage to perform conversion
  allocate(tmp_array(nglob),stat=ier)
  if( ier /=0 ) call exit_MPI(myrank,'error allocating temporary array for mesh topology')

  !--- x coordinate
  tmp_array(:) = 0._CUSTOM_REAL
  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          ! distinguish between single and double precision for reals
          if(CUSTOM_REAL == SIZE_REAL) then
            tmp_array(iglob) = sngl(xstore(i,j,k,ispec))
          else
            tmp_array(iglob) = xstore(i,j,k,ispec)
          endif
        enddo
      enddo
    enddo
  enddo
  write(27) tmp_array ! xstore

  !--- y coordinate
  tmp_array(:) = 0._CUSTOM_REAL
  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          ! distinguish between single and double precision for reals
          if(CUSTOM_REAL == SIZE_REAL) then
            tmp_array(iglob) = sngl(ystore(i,j,k,ispec))
          else
            tmp_array(iglob) = ystore(i,j,k,ispec)
          endif
        enddo
      enddo
    enddo
  enddo
  write(27) tmp_array ! ystore

  !--- z coordinate
  tmp_array(:) = 0._CUSTOM_REAL
  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          ! distinguish between single and double precision for reals
          if(CUSTOM_REAL == SIZE_REAL) then
            tmp_array(iglob) = sngl(zstore(i,j,k,ispec))
          else
            tmp_array(iglob) = zstore(i,j,k,ispec)
          endif
        enddo
      enddo
    enddo
  enddo
  write(27) tmp_array ! zstore

  deallocate(tmp_array)

  ! local to global indexing
  write(27) ibool
  write(27) idoubling
  write(27) ispec_is_tiso

  ! local GLL points
  write(27) xixstore
  write(27) xiystore
  write(27) xizstore
  write(27) etaxstore
  write(27) etaystore
  write(27) etazstore
  write(27) gammaxstore
  write(27) gammaystore
  write(27) gammazstore

  write(27) rhostore
  write(27) kappavstore

  ! other terms needed in the solid regions only
  if(iregion_code /= IREGION_OUTER_CORE) then

    ! note: muvstore needed for Q_mu shear attenuation in inner core
    if(.not. (ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE)) then
      write(27) muvstore
    endif

    !   save anisotropy in the mantle only
    if(TRANSVERSE_ISOTROPY) then
      if(iregion_code == IREGION_CRUST_MANTLE .and. .not. ANISOTROPIC_3D_MANTLE) then
        write(27) kappahstore
        write(27) muhstore
        write(27) eta_anisostore
      endif
    endif

    !   save anisotropy in the inner core only
    if(ANISOTROPIC_INNER_CORE .and. iregion_code == IREGION_INNER_CORE) then
      write(27) c11store
      write(27) c33store
      write(27) c12store
      write(27) c13store
      write(27) c44store
    endif

    if(ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
        write(27) c11store
        write(27) c12store
        write(27) c13store
        write(27) c14store
        write(27) c15store
        write(27) c16store
        write(27) c22store
        write(27) c23store
        write(27) c24store
        write(27) c25store
        write(27) c26store
        write(27) c33store
        write(27) c34store
        write(27) c35store
        write(27) c36store
        write(27) c44store
        write(27) c45store
        write(27) c46store
        write(27) c55store
        write(27) c56store
        write(27) c66store
    endif

  endif

  ! Stacey
  if(ABSORBING_CONDITIONS) then

    if(iregion_code == IREGION_CRUST_MANTLE) then
      write(27) rho_vp
      write(27) rho_vs
    else if(iregion_code == IREGION_OUTER_CORE) then
      write(27) rho_vp
    endif

  endif

  ! mass matrices
  !
  ! in the case of stacey boundary conditions, add C*deltat/2 contribution to the mass matrix
  ! on the Stacey edges for the crust_mantle and outer_core regions but not for the inner_core region
  ! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix
  !
  ! if absorbing_conditions are not set or if NCHUNKS=6, only one mass matrix is needed
  ! for the sake of performance, only "rmassz" array will be filled and "rmassx" & "rmassy" will be obsolete
  if(NCHUNKS /= 6 .and. ABSORBING_CONDITIONS .and. iregion_code == IREGION_CRUST_MANTLE) then
     write(27) rmassx
     write(27) rmassy
  endif

  write(27) rmassz

  ! additional ocean load mass matrix if oceans and if we are in the crust
  if(OCEANS .and. iregion_code == IREGION_CRUST_MANTLE) write(27) rmass_ocean_load

  close(27) ! solver_data.bin

  ! absorbing boundary parameters
  open(unit=27,file=prname(1:len_trim(prname))//'boundary.bin', &
        status='unknown',form='unformatted',action='write',iostat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error opening boundary.bin file')

  write(27) nspec2D_xmin
  write(27) nspec2D_xmax
  write(27) nspec2D_ymin
  write(27) nspec2D_ymax
  write(27) NSPEC2D_BOTTOM
  write(27) NSPEC2D_TOP

  write(27) ibelm_xmin
  write(27) ibelm_xmax
  write(27) ibelm_ymin
  write(27) ibelm_ymax
  write(27) ibelm_bottom
  write(27) ibelm_top

  write(27) normal_xmin
  write(27) normal_xmax
  write(27) normal_ymin
  write(27) normal_ymax
  write(27) normal_bottom
  write(27) normal_top

  write(27) jacobian2D_xmin
  write(27) jacobian2D_xmax
  write(27) jacobian2D_ymin
  write(27) jacobian2D_ymax
  write(27) jacobian2D_bottom
  write(27) jacobian2D_top

  close(27)

  if(ATTENUATION) then
    open(unit=27, file=prname(1:len_trim(prname))//'attenuation.bin', &
          status='unknown', form='unformatted',action='write',iostat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error opening attenuation.bin file')

    write(27) tau_s
    write(27) tau_e_store
    write(27) Qmu_store
    write(27) T_c_source
    close(27)
  endif

  if(HETEROGEN_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
    open(unit=27,file=prname(1:len_trim(prname))//'dvp.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error opening dvp.bin file')

    write(27) dvpstore
    close(27)
  endif

  ! uncomment for vp & vs model storage
  if( SAVE_MESH_FILES ) then
    ! outputs model files in binary format
    call save_arrays_solver_meshfiles(myrank,nspec)
  endif

  end subroutine save_arrays_solver

!
!-------------------------------------------------------------------------------------------------
!

  subroutine save_arrays_solver_meshfiles(myrank,nspec)

! outputs model files in binary format

  use constants

  use meshfem3D_models_par,only: &
    TRANSVERSE_ISOTROPY,ATTENUATION

  use create_regions_mesh_par2,only: &
    rhostore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
    Qmu_store, &
    prname

  implicit none

  integer :: myrank
  integer :: nspec

  ! local parameters
  integer :: i,j,k,ispec,ier
  real(kind=CUSTOM_REAL) :: scaleval1,scaleval2
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: temp_store

  ! scaling factors to re-dimensionalize units
  scaleval1 = sngl( sqrt(PI*GRAV*RHOAV)*(R_EARTH/1000.0d0) )
  scaleval2 = sngl( RHOAV/1000.0d0 )

  ! isotropic model
  ! vp
  open(unit=27,file=prname(1:len_trim(prname))//'vp.bin', &
       status='unknown',form='unformatted',action='write',iostat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error opening vp.bin file')

  write(27) sqrt( (kappavstore+4.*muvstore/3.)/rhostore )*scaleval1
  close(27)
  ! vs
  open(unit=27,file=prname(1:len_trim(prname))//'vs.bin', &
        status='unknown',form='unformatted',action='write',iostat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error opening vs.bin file')

  write(27) sqrt( muvstore/rhostore )*scaleval1
  close(27)
  ! rho
  open(unit=27,file=prname(1:len_trim(prname))//'rho.bin', &
        status='unknown',form='unformatted',action='write',iostat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error opening rho.bin file')

  write(27) rhostore*scaleval2
  close(27)

  ! transverse isotropic model
  if( TRANSVERSE_ISOTROPY ) then
    ! vpv
    open(unit=27,file=prname(1:len_trim(prname))//'vpv.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error opening vpv.bin file')

    write(27) sqrt( (kappavstore+4.*muvstore/3.)/rhostore )*scaleval1
    close(27)
    ! vph
    open(unit=27,file=prname(1:len_trim(prname))//'vph.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error opening vph.bin file')

    write(27) sqrt( (kappahstore+4.*muhstore/3.)/rhostore )*scaleval1
    close(27)
    ! vsv
    open(unit=27,file=prname(1:len_trim(prname))//'vsv.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error opening vsv.bin file')

    write(27) sqrt( muvstore/rhostore )*scaleval1
    close(27)
    ! vsh
    open(unit=27,file=prname(1:len_trim(prname))//'vsh.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error opening vsh.bin file')

    write(27) sqrt( muhstore/rhostore )*scaleval1
    close(27)
    ! rho
    open(unit=27,file=prname(1:len_trim(prname))//'rho.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error opening rho.bin file')

    write(27) rhostore*scaleval2
    close(27)
    ! eta
    open(unit=27,file=prname(1:len_trim(prname))//'eta.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error opening eta.bin file')

    write(27) eta_anisostore
    close(27)
  endif ! TRANSVERSE_ISOTROPY

  ! shear attenuation
  if( ATTENUATION ) then
    ! saves Qmu_store to full custom_real array
    ! uses temporary array
    allocate(temp_store(NGLLX,NGLLY,NGLLZ,nspec))
    if (USE_3D_ATTENUATION_ARRAYS) then
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

    ! Qmu
    open(unit=27,file=prname(1:len_trim(prname))//'qmu.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error opening qmu.bin file')

    write(27) temp_store
    close(27)

    ! frees temporary memory
    deallocate(temp_store)
  endif ! ATTENUATION

  end subroutine save_arrays_solver_meshfiles
!
!-------------------------------------------------------------------------------------------------
!


  subroutine save_arrays_solver_MPI(iregion_code)

  use meshfem3D_par,only: &
    myrank,LOCAL_PATH, &
    IREGION_CRUST_MANTLE,IREGION_OUTER_CORE,IREGION_INNER_CORE, &
    ADIOS_FOR_MPI_ARRAYS

!  use create_MPI_interfaces_par

  use MPI_crust_mantle_par
  use MPI_outer_core_par
  use MPI_inner_core_par

  implicit none

  integer,intent(in):: iregion_code

  select case( iregion_code )
  case( IREGION_CRUST_MANTLE )
    ! crust mantle
    if (ADIOS_FOR_MPI_ARRAYS) then
      call save_MPI_arrays_adios(myrank,IREGION_CRUST_MANTLE,LOCAL_PATH, &
          num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
          my_neighbours_crust_mantle,nibool_interfaces_crust_mantle, &
          ibool_interfaces_crust_mantle, &
          nspec_inner_crust_mantle,nspec_outer_crust_mantle, &
          num_phase_ispec_crust_mantle,phase_ispec_inner_crust_mantle, &
          num_colors_outer_crust_mantle,num_colors_inner_crust_mantle, &
          num_elem_colors_crust_mantle)
    else
      call save_MPI_arrays(myrank,IREGION_CRUST_MANTLE,LOCAL_PATH, &
          num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
          my_neighbours_crust_mantle,nibool_interfaces_crust_mantle, &
          ibool_interfaces_crust_mantle, &
          nspec_inner_crust_mantle,nspec_outer_crust_mantle, &
          num_phase_ispec_crust_mantle,phase_ispec_inner_crust_mantle, &
          num_colors_outer_crust_mantle,num_colors_inner_crust_mantle, &
          num_elem_colors_crust_mantle)
    endif

  case( IREGION_OUTER_CORE )
    ! outer core
    if (ADIOS_FOR_MPI_ARRAYS) then
      call save_MPI_arrays_adios(myrank,IREGION_OUTER_CORE,LOCAL_PATH, &
          num_interfaces_outer_core,max_nibool_interfaces_oc, &
          my_neighbours_outer_core,nibool_interfaces_outer_core, &
          ibool_interfaces_outer_core, &
          nspec_inner_outer_core,nspec_outer_outer_core, &
          num_phase_ispec_outer_core,phase_ispec_inner_outer_core, &
          num_colors_outer_outer_core,num_colors_inner_outer_core, &
          num_elem_colors_outer_core)
    else
      call save_MPI_arrays(myrank,IREGION_OUTER_CORE,LOCAL_PATH, &
          num_interfaces_outer_core,max_nibool_interfaces_oc, &
          my_neighbours_outer_core,nibool_interfaces_outer_core, &
          ibool_interfaces_outer_core, &
          nspec_inner_outer_core,nspec_outer_outer_core, &
          num_phase_ispec_outer_core,phase_ispec_inner_outer_core, &
          num_colors_outer_outer_core,num_colors_inner_outer_core, &
          num_elem_colors_outer_core)
    endif

  case( IREGION_INNER_CORE )
    ! inner core
    if (ADIOS_FOR_MPI_ARRAYS) then
      call save_MPI_arrays_adios(myrank,IREGION_INNER_CORE,LOCAL_PATH, &
          num_interfaces_inner_core,max_nibool_interfaces_ic, &
          my_neighbours_inner_core,nibool_interfaces_inner_core, &
          ibool_interfaces_inner_core, &
          nspec_inner_inner_core,nspec_outer_inner_core, &
          num_phase_ispec_inner_core,phase_ispec_inner_inner_core, &
          num_colors_outer_inner_core,num_colors_inner_inner_core, &
          num_elem_colors_inner_core)
    else
      call save_MPI_arrays(myrank,IREGION_INNER_CORE,LOCAL_PATH, &
          num_interfaces_inner_core,max_nibool_interfaces_ic, &
          my_neighbours_inner_core,nibool_interfaces_inner_core, &
          ibool_interfaces_inner_core, &
          nspec_inner_inner_core,nspec_outer_inner_core, &
          num_phase_ispec_inner_core,phase_ispec_inner_inner_core, &
          num_colors_outer_inner_core,num_colors_inner_inner_core, &
          num_elem_colors_inner_core)
    endif
  end select

  end subroutine save_arrays_solver_MPI


!
!-------------------------------------------------------------------------------------------------
!

  subroutine save_MPI_arrays(myrank,iregion_code,LOCAL_PATH, &
                                  num_interfaces,max_nibool_interfaces, &
                                  my_neighbours,nibool_interfaces, &
                                  ibool_interfaces, &
                                  nspec_inner,nspec_outer, &
                                  num_phase_ispec,phase_ispec_inner, &
                                  num_colors_outer,num_colors_inner, &
                                  num_elem_colors)
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
  character(len=150) :: prname
  integer :: ier

  ! create the name for the database of the current slide and region
  call create_name_database(prname,myrank,iregion_code,LOCAL_PATH)

  open(unit=IOUT,file=prname(1:len_trim(prname))//'solver_data_mpi.bin', &
       status='unknown',action='write',form='unformatted',iostat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error opening solver_data_mpi.bin')

  ! MPI interfaces
  write(IOUT) num_interfaces
  if( num_interfaces > 0 ) then
    write(IOUT) max_nibool_interfaces
    write(IOUT) my_neighbours
    write(IOUT) nibool_interfaces
    write(IOUT) ibool_interfaces
  endif

  ! inner/outer elements
  write(IOUT) nspec_inner,nspec_outer
  write(IOUT) num_phase_ispec
  if(num_phase_ispec > 0 ) write(IOUT) phase_ispec_inner

  ! mesh coloring
  if( USE_MESH_COLORING_GPU ) then
    write(IOUT) num_colors_outer,num_colors_inner
    write(IOUT) num_elem_colors
  endif

  close(IOUT)

  end subroutine save_MPI_arrays


!
!-------------------------------------------------------------------------------------------------
!

  subroutine save_arrays_solver_boundary()

! saves arrays for boundaries such as MOHO, 400 and 670 discontinuities

  use meshfem3d_par,only: &
    myrank

  use meshfem3D_models_par,only: &
    SAVE_BOUNDARY_MESH,HONOR_1D_SPHERICAL_MOHO,SUPPRESS_CRUSTAL_MESH

  use create_regions_mesh_par2, only: &
    NSPEC2D_MOHO, NSPEC2D_400, NSPEC2D_670, &
    ibelm_moho_top,ibelm_moho_bot,ibelm_400_top,ibelm_400_bot, &
    ibelm_670_top,ibelm_670_bot,normal_moho,normal_400,normal_670, &
    ispec2D_moho_top,ispec2D_moho_bot,ispec2D_400_top,ispec2D_400_bot, &
    ispec2D_670_top,ispec2D_670_bot, &
    prname

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

  end subroutine save_arrays_solver_boundary

