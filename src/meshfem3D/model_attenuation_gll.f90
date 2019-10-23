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

!--------------------------------------------------------------------------------------------------
! GLL QMU
!
! based on modified GLL mesh output from mesher
!
! used for iterative inversion procedures
!--------------------------------------------------------------------------------------------------

  module model_gll_qmu_par

  use constants, only: CUSTOM_REAL

  ! GLL model_variables
  type model_gll_qmu_variables
    sequence
    ! tomographic iteration model on GLL points

    real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: qmu_new

    ! number of elements (crust/mantle elements)
    integer :: nspec
    integer :: dummy_pad ! padding 4 bytes to align the structure
  end type model_gll_qmu_variables
  type (model_gll_qmu_variables) :: MGLL_V

  ! model GLL type: 1 == iso, 2 == tiso, 3 == azi
  integer :: MGLL_TYPE

  end module model_gll_qmu_par

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_attenuation_gll_broadcast()

! standard routine to setup model

  use constants

  use shared_parameters, only: NSPEC_REGIONS,ADIOS_FOR_MODELS,NPROCTOT,NCHUNKS, &
                               MODEL,MODEL_GLL_TYPE

  use model_gll_qmu_par

  implicit none

  ! local parameters
  double precision :: scaleval
  real(kind=CUSTOM_REAL) :: minvalue,maxvalue,min_all,max_all
  integer :: ier

  ! sets type (iso,tiso,azi)
  MGLL_TYPE = MODEL_GLL_TYPE

  ! sets number of elements (crust/mantle region)
  MGLL_V%nspec = NSPEC_REGIONS(IREGION_CRUST_MANTLE)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'broadcast model: ',trim(MODEL)
    write(IMAIN,*) 'setup:'
    write(IMAIN,*) '  NCHUNKS           : ',NCHUNKS
    write(IMAIN,*) '  NPROC total       : ',NPROCTOT
    write(IMAIN,*) '  NSPEC             : ',MGLL_V%nspec
    write(IMAIN,*) '  NGLLX/NGLLY/NGLLZ : ',NGLLX,NGLLY,NGLLZ
    call flush_IMAIN()
  endif

  ! safety check
  if (MGLL_TYPE < 1 .or. MGLL_TYPE > 3) &
    stop 'Invalid MODEL_GLL_TYPE, please use 1(iso), 2(tiso) or 3(azi) in get_model_parameters.F90 setup'

  ! allocates arrays
  ! differs for isotropic model or transverse isotropic models
  allocate( MGLL_V%qmu_new(NGLLX,NGLLY,NGLLZ,MGLL_V%nspec), &
       stat=ier)
  MGLL_V%qmu_new(:,:,:,:) = 0.0_CUSTOM_REAL
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating qmu_new,.. arrays')

  ! reads in model files for each process
  if (ADIOS_FOR_MODELS) then
     call read_gll_qmu_model_adios(myrank)
  else
     call exit_MPI(myrank,'GLL Qmu models are not implemented for binary files.')
  endif

  call print_min_max_all(MGLL_V%qmu_new,"qmu new")

  call synchronize_all()

  contains

    subroutine print_min_max_all(array,name)

    implicit none

    real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,MGLL_V%nspec),intent(in) :: array
    character(len=*),intent(in) :: name
    ! local parameters
    real(kind=CUSTOM_REAL) :: minvalue,maxvalue,min_all,max_all

    maxvalue = maxval( array )
    minvalue = minval( array )
    call max_all_cr(maxvalue, max_all)
    call min_all_cr(minvalue, min_all)

    if (myrank == 0) then
      write(IMAIN,*) '  '//trim(name)//' min/max: ',min_all,max_all
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    end subroutine

  end subroutine model_attenuation_gll_broadcast


  
  subroutine model_attenuation_gll(ispec, i, j, k, Qmu)

  use model_gll_qmu_par

  implicit none

  integer,intent(in) :: ispec,i,j,k
  double precision,intent(inout) :: Qmu

  Qmu = dble( MGLL_V%qmu_new(i,j,k,ispec) )

  end subroutine model_attenuation_gll
