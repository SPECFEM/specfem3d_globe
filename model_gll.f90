!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            March 2010
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

!--------------------------------------------------------------------------------------------------
! GLL
!
! based on modified GLL mesh output from mesher
!
! used for iterative inversion procedures
!--------------------------------------------------------------------------------------------------


  subroutine model_gll_broadcast(myrank,MGLL_V,NSPEC)

! standard routine to setup model

  use meshfem3D_models_par,only: TRANSVERSE_ISOTROPY

  implicit none

  include "constants.h"
  ! standard include of the MPI library
  include 'mpif.h'

  ! GLL model_variables
  type model_gll_variables
    sequence
    ! tomographic iteration model on GLL points
    double precision :: scale_velocity,scale_density
    ! isotropic model
    real(kind=CUSTOM_REAL),dimension(:,:,:,:),pointer :: vs_new,vp_new,rho_new
    ! transverse isotropic model
    real(kind=CUSTOM_REAL),dimension(:,:,:,:),pointer :: vsv_new,vpv_new, &
      vsh_new,vph_new,eta_new
    logical :: MODEL_GLL
    logical,dimension(3) :: dummy_pad ! padding 3 bytes to align the structure
  end type model_gll_variables
  type (model_gll_variables) MGLL_V

  integer, dimension(MAX_NUM_REGIONS) :: NSPEC
  integer :: myrank

  ! local parameters
  double precision :: min_dvs,max_dvs,min_dvs_all,max_dvs_all,scaleval
  integer :: ier

  ! differs for isotropic model or transverse isotropic models
  if( .not. TRANSVERSE_ISOTROPY ) then 
    ! isotropic model
    allocate( MGLL_V%vp_new(NGLLX,NGLLY,NGLLZ,NSPEC(IREGION_CRUST_MANTLE)) )
    allocate( MGLL_V%vs_new(NGLLX,NGLLY,NGLLZ,NSPEC(IREGION_CRUST_MANTLE)) )
  else
    ! transverse isotropic model
    allocate( MGLL_V%vpv_new(NGLLX,NGLLY,NGLLZ,NSPEC(IREGION_CRUST_MANTLE)) )
    allocate( MGLL_V%vsv_new(NGLLX,NGLLY,NGLLZ,NSPEC(IREGION_CRUST_MANTLE)) )
    allocate( MGLL_V%vph_new(NGLLX,NGLLY,NGLLZ,NSPEC(IREGION_CRUST_MANTLE)) )
    allocate( MGLL_V%vsh_new(NGLLX,NGLLY,NGLLZ,NSPEC(IREGION_CRUST_MANTLE)) )
    allocate( MGLL_V%eta_new(NGLLX,NGLLY,NGLLZ,NSPEC(IREGION_CRUST_MANTLE)) )        
  endif
  allocate( MGLL_V%rho_new(NGLLX,NGLLY,NGLLZ,NSPEC(IREGION_CRUST_MANTLE)) )
  
  ! non-dimensionalize scaling values
  ! (model velocities must be given as km/s)
  scaleval = dsqrt(PI*GRAV*RHOAV)
  MGLL_V%scale_velocity = 1000.0d0/(R_EARTH*scaleval)
  MGLL_V%scale_density =  1000.0d0/RHOAV

  call read_gll_model(myrank,MGLL_V,NSPEC)

  ! checks velocity range
  if( .not. TRANSVERSE_ISOTROPY ) then       
    max_dvs = maxval( MGLL_V%vs_new )
    min_dvs = minval( MGLL_V%vs_new )
  else
    max_dvs = maxval( MGLL_V%vsv_new + MGLL_V%vsh_new )
    min_dvs = minval( MGLL_V%vsv_new + MGLL_V%vsh_new )  
  endif
  call mpi_reduce(max_dvs, max_dvs_all, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD,ier)
  call mpi_reduce(min_dvs, min_dvs_all, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD,ier)
  if( myrank == 0 ) then
    write(IMAIN,*)'model GLL:'
    write(IMAIN,*) '  vs new min/max: ',min_dvs_all,max_dvs_all
    write(IMAIN,*)
  endif

  end subroutine model_gll_broadcast

!
!-------------------------------------------------------------------------------------------------
!


  subroutine read_gll_model(myrank,MGLL_V,NSPEC)

  use meshfem3D_models_par,only: TRANSVERSE_ISOTROPY
  
  implicit none

  include "constants.h"

  ! GLL model_variables
  type model_gll_variables
    sequence
    ! tomographic iteration model on GLL points
    double precision :: scale_velocity,scale_density
    ! isotropic model
    real(kind=CUSTOM_REAL),dimension(:,:,:,:),pointer :: vs_new,vp_new,rho_new
    ! transverse isotropic model
    real(kind=CUSTOM_REAL),dimension(:,:,:,:),pointer :: vsv_new,vpv_new, &
      vsh_new,vph_new,eta_new
    logical :: MODEL_GLL
    logical,dimension(3) :: dummy_pad ! padding 3 bytes to align the structure
  end type model_gll_variables
  type (model_gll_variables) MGLL_V

  integer, dimension(MAX_NUM_REGIONS) :: NSPEC
  integer :: myrank

  !--------------------------------------------------------------------
  ! USER PARAMETER

  character(len=150),parameter:: MGLL_path = 'DATA/GLL/'
  !--------------------------------------------------------------------

  ! local parameters
  integer :: ier
  character(len=150) :: prname

  if( myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*)'reading in model from ',trim(MGLL_path)
  endif
  
  ! only crust and mantle
  write(prname,'(a,i6.6,a)') MGLL_path(1:len_trim(MGLL_path))//'proc',myrank,'_reg1_'

  ! reads in model for each partition
  if( .not. TRANSVERSE_ISOTROPY ) then       
    ! isotropic model
    ! vp mesh
    open(unit=27,file=prname(1:len_trim(prname))//'vp_new.bin',&
          status='old',action='read',form='unformatted',iostat=ier)
    if( ier /= 0 ) then
      write(IMAIN,*) 'error opening: ',prname(1:len_trim(prname))//'vp_new.bin'
      call exit_MPI(myrank,'error model gll')
    endif
    read(27) MGLL_V%vp_new(:,:,:,1:nspec(IREGION_CRUST_MANTLE))
    close(27)

    ! vs mesh
    open(unit=27,file=prname(1:len_trim(prname))//'vs_new.bin', &
         status='old',action='read',form='unformatted',iostat=ier)
    if( ier /= 0 ) then
      print*,'error opening: ',prname(1:len_trim(prname))//'vs_new.bin'
      call exit_MPI(myrank,'error model gll')
    endif
    read(27) MGLL_V%vs_new(:,:,:,1:nspec(IREGION_CRUST_MANTLE))
    close(27)

  else
  
    ! transverse isotropic model
    ! vp mesh
    open(unit=27,file=prname(1:len_trim(prname))//'vpv_new.bin',&
          status='old',action='read',form='unformatted',iostat=ier)
    if( ier /= 0 ) then
      write(IMAIN,*) 'error opening: ',prname(1:len_trim(prname))//'vpv_new.bin'
      call exit_MPI(myrank,'error model gll')
    endif
    read(27) MGLL_V%vpv_new(:,:,:,1:nspec(IREGION_CRUST_MANTLE))
    close(27)

    open(unit=27,file=prname(1:len_trim(prname))//'vph_new.bin',&
          status='old',action='read',form='unformatted',iostat=ier)
    if( ier /= 0 ) then
      write(IMAIN,*) 'error opening: ',prname(1:len_trim(prname))//'vph_new.bin'
      call exit_MPI(myrank,'error model gll')
    endif
    read(27) MGLL_V%vph_new(:,:,:,1:nspec(IREGION_CRUST_MANTLE))
    close(27)

    ! vs mesh
    open(unit=27,file=prname(1:len_trim(prname))//'vsv_new.bin', &
         status='old',action='read',form='unformatted',iostat=ier)
    if( ier /= 0 ) then
      print*,'error opening: ',prname(1:len_trim(prname))//'vsv_new.bin'
      call exit_MPI(myrank,'error model gll')
    endif
    read(27) MGLL_V%vsv_new(:,:,:,1:nspec(IREGION_CRUST_MANTLE))
    close(27)

    open(unit=27,file=prname(1:len_trim(prname))//'vsh_new.bin', &
         status='old',action='read',form='unformatted',iostat=ier)
    if( ier /= 0 ) then
      print*,'error opening: ',prname(1:len_trim(prname))//'vsh_new.bin'
      call exit_MPI(myrank,'error model gll')
    endif
    read(27) MGLL_V%vsh_new(:,:,:,1:nspec(IREGION_CRUST_MANTLE))
    close(27)

    ! eta mesh
    open(unit=27,file=prname(1:len_trim(prname))//'eta_new.bin', &
         status='old',action='read',form='unformatted',iostat=ier)
    if( ier /= 0 ) then
      print*,'error opening: ',prname(1:len_trim(prname))//'eta_new.bin'
      call exit_MPI(myrank,'error model gll')
    endif
    read(27) MGLL_V%eta_new(:,:,:,1:nspec(IREGION_CRUST_MANTLE))
    close(27)

  endif

  ! rho mesh
  open(unit=27,file=prname(1:len_trim(prname))//'rho_new.bin', &
       status='old',action='read',form='unformatted',iostat=ier)
  if( ier /= 0 ) then
    print*,'error opening: ',prname(1:len_trim(prname))//'rho_new.bin'
    call exit_MPI(myrank,'error model gll')
  endif
  read(27) MGLL_V%rho_new(:,:,:,1:nspec(IREGION_CRUST_MANTLE))
  close(27)

  end subroutine read_gll_model
