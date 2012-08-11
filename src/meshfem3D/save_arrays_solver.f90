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

  subroutine save_arrays_solver(myrank,nspec,nglob,idoubling,ibool, &
                    iregion_code,xstore,ystore,zstore, &
                    is_on_a_slice_edge, &
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

  ! this for non blocking MPI
  logical, dimension(nspec) :: is_on_a_slice_edge

  ! boundary parameters locator
  integer :: NSPEC2D_TOP,NSPEC2D_BOTTOM

  ! local parameters
  integer i,j,k,ispec,iglob,nspec1,nglob1,ier
  real(kind=CUSTOM_REAL) :: scaleval1,scaleval2

  ! save nspec and nglob, to be used in combine_paraview_data
  open(unit=27,file=prname(1:len_trim(prname))//'array_dims.txt', &
        status='unknown',action='write',iostat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error opening array_dims file')

  nspec1 = nspec
  nglob1 = nglob

  write(27,*) nspec1
  write(27,*) nglob1
  close(27)

  ! mesh parameters
  open(unit=27,file=prname(1:len_trim(prname))//'solver_data_1.bin', &
       status='unknown',form='unformatted',action='write',iostat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error opening solver_data_1.bin file')

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

  if(HETEROGEN_3D_MANTLE) then
     open(unit=29,file=prname(1:len_trim(prname))//'dvp.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
     if( ier /= 0 ) call exit_mpi(myrank,'error opening dvp.bin file')

     write(29) dvpstore
     close(29)
  endif

  ! other terms needed in the solid regions only
  if(iregion_code /= IREGION_OUTER_CORE) then

    if(.not. (ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE)) write(27) muvstore

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
  ! on Stacey edges for the crust_mantle and outer_core regions but not for the inner_core region
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

  close(27)

  ! mesh topology
  open(unit=27,file=prname(1:len_trim(prname))//'solver_data_2.bin', &
        status='unknown',form='unformatted',action='write',iostat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error opening solver_data_2.bin file')

! mesh arrays used in the solver to locate source and receivers
! and for anisotropy and gravity, save in single precision
! use rmassz for temporary storage to perform conversion, since already saved

  !--- x coordinate
  rmassz(:) = 0._CUSTOM_REAL
  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          ! distinguish between single and double precision for reals
          if(CUSTOM_REAL == SIZE_REAL) then
            rmassz(iglob) = sngl(xstore(i,j,k,ispec))
          else
            rmassz(iglob) = xstore(i,j,k,ispec)
          endif
        enddo
      enddo
    enddo
  enddo
  write(27) rmassz

  !--- y coordinate
  rmassz(:) = 0._CUSTOM_REAL
  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          ! distinguish between single and double precision for reals
          if(CUSTOM_REAL == SIZE_REAL) then
            rmassz(iglob) = sngl(ystore(i,j,k,ispec))
          else
            rmassz(iglob) = ystore(i,j,k,ispec)
          endif
        enddo
      enddo
    enddo
  enddo
  write(27) rmassz

  !--- z coordinate
  rmassz(:) = 0._CUSTOM_REAL
  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          ! distinguish between single and double precision for reals
          if(CUSTOM_REAL == SIZE_REAL) then
            rmassz(iglob) = sngl(zstore(i,j,k,ispec))
          else
            rmassz(iglob) = zstore(i,j,k,ispec)
          endif
        enddo
      enddo
    enddo
  enddo
  write(27) rmassz

  write(27) ibool
  write(27) idoubling
  write(27) is_on_a_slice_edge
  write(27) ispec_is_tiso

  close(27)

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

  ! uncomment for vp & vs model storage
  if( SAVE_MESH_FILES ) then
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
    endif
  endif ! SAVE_MESH_FILES

  end subroutine save_arrays_solver

!
!-------------------------------------------------------------------------------------------------
!


  subroutine save_arrays_solver_MPI(iregion_code)

  use meshfem3D_par,only: &
    myrank,LOCAL_PATH, &
    IREGION_CRUST_MANTLE,IREGION_OUTER_CORE,IREGION_INNER_CORE

!  use create_MPI_interfaces_par

  use MPI_crust_mantle_par
  use MPI_outer_core_par
  use MPI_inner_core_par

  implicit none

  integer,intent(in):: iregion_code

  select case( iregion_code )
  case( IREGION_CRUST_MANTLE )
    ! crust mantle
    call save_MPI_arrays(myrank,IREGION_CRUST_MANTLE,LOCAL_PATH, &
                             num_interfaces_crust_mantle,max_nibool_interfaces_crust_mantle, &
                             my_neighbours_crust_mantle,nibool_interfaces_crust_mantle, &
                             ibool_interfaces_crust_mantle, &
                             nspec_inner_crust_mantle,nspec_outer_crust_mantle, &
                             num_phase_ispec_crust_mantle,phase_ispec_inner_crust_mantle, &
                             num_colors_outer_crust_mantle,num_colors_inner_crust_mantle, &
                             num_elem_colors_crust_mantle)


  case( IREGION_OUTER_CORE )
    ! outer core
    call save_MPI_arrays(myrank,IREGION_OUTER_CORE,LOCAL_PATH, &
                             num_interfaces_outer_core,max_nibool_interfaces_outer_core, &
                             my_neighbours_outer_core,nibool_interfaces_outer_core, &
                             ibool_interfaces_outer_core, &
                             nspec_inner_outer_core,nspec_outer_outer_core, &
                             num_phase_ispec_outer_core,phase_ispec_inner_outer_core, &
                             num_colors_outer_outer_core,num_colors_inner_outer_core, &
                             num_elem_colors_outer_core)

  case( IREGION_INNER_CORE )
    ! inner core
    call save_MPI_arrays(myrank,IREGION_INNER_CORE,LOCAL_PATH, &
                             num_interfaces_inner_core,max_nibool_interfaces_inner_core, &
                             my_neighbours_inner_core,nibool_interfaces_inner_core, &
                             ibool_interfaces_inner_core, &
                             nspec_inner_inner_core,nspec_outer_inner_core, &
                             num_phase_ispec_inner_core,phase_ispec_inner_inner_core, &
                             num_colors_outer_inner_core,num_colors_inner_inner_core, &
                             num_elem_colors_inner_core)

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

