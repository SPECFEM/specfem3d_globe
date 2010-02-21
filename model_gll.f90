!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
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
  
  implicit none

  include "constants.h"
  ! standard include of the MPI library
  include 'mpif.h'

  ! GLL model_variables
  type model_gll_variables
    ! tomographic iteration model on GLL points
    real(kind=CUSTOM_REAL),dimension(:,:,:,:),pointer :: vs_new,vp_new,rho_new
    double precision :: scale_velocity,scale_density
    logical :: MODEL_GLL  
  end type model_gll_variables
  type (model_gll_variables) MGLL_V

  integer, dimension(MAX_NUM_REGIONS) :: NSPEC
  integer :: myrank
    
  ! local parameters
  double precision :: min_dvs,max_dvs
  integer :: ier

  allocate( MGLL_V%vp_new(NGLLX,NGLLY,NGLLZ,NSPEC(IREGION_CRUST_MANTLE)) )
  allocate( MGLL_V%vs_new(NGLLX,NGLLY,NGLLZ,NSPEC(IREGION_CRUST_MANTLE)) )
  allocate( MGLL_V%rho_new(NGLLX,NGLLY,NGLLZ,NSPEC(IREGION_CRUST_MANTLE)) )
  ! non-dimensionalize scaling values
  MGLL_V%scale_velocity = 1000.0d0/(PI*GRAV*RHOAV*R_EARTH)
  MGLL_V%scale_density =  1000.0d0/RHOAV
  
  call read_gll_model(myrank,MGLL_V,NSPEC)
  
  ! checks velocity range 
  max_dvs = maxval( MGLL_V%vs_new )
  min_dvs = minval( MGLL_V%vs_new )
  call mpi_reduce(max_dvs, max_dvs, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD,ier)
  call mpi_reduce(min_dvs, min_dvs, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD,ier)
  if( myrank == 0 ) then
    write(IMAIN,*)'model GLL:'
    write(IMAIN,*) '  vs new min/max: ',min_dvs,max_dvs    
    write(IMAIN,*)
  endif      
  
  end subroutine model_gll_broadcast

!
!-------------------------------------------------------------------------------------------------
!
  
  
  subroutine read_gll_model(myrank,MGLL_V,NSPEC)

  implicit none

  include "constants.h"
  
  ! GLL model_variables
  type model_gll_variables
    ! tomographic iteration model on GLL points
    real(kind=CUSTOM_REAL),dimension(:,:,:,:),pointer :: vs_new,vp_new,rho_new
    double precision :: scale_velocity,scale_density
    logical :: MODEL_GLL  
  end type model_gll_variables
  type (model_gll_variables) MGLL_V  

  integer, dimension(MAX_NUM_REGIONS) :: NSPEC  
  integer :: myrank

  !--------------------------------------------------------------------
  ! USER PARAMETER
  
  character(len=150),parameter:: MGLL_path = 'KERNELS/model_m1/'      
  !--------------------------------------------------------------------

  ! local parameters
  integer :: ier
  character(len=150) :: prname
  
  ! only crust and mantle
  write(prname,'(a,i6.6,a)') MGLL_path(1:len_trim(MGLL_path))//'proc',myrank,'_reg1_'
      
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