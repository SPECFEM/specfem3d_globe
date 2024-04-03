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
  type (model_gll_qmu_variables) :: MGLL_QMU_V

  end module model_gll_qmu_par

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_attenuation_gll_broadcast()

! standard routine to setup model

  use constants

  use shared_parameters, only: NSPEC_REGIONS,ADIOS_FOR_MODELS,NPROCTOT,NCHUNKS, &
                               MODEL

  use model_gll_par, only: MGLL_V
  use model_gll_qmu_par

  implicit none

  ! local parameters
  integer :: ier

  ! sets number of elements (crust/mantle region)
  MGLL_QMU_V%nspec = NSPEC_REGIONS(IREGION_CRUST_MANTLE)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'broadcast model: ',trim(MODEL)
    write(IMAIN,*) 'setup:'
    write(IMAIN,*) '  NCHUNKS           : ',NCHUNKS
    write(IMAIN,*) '  NPROC total       : ',NPROCTOT
    write(IMAIN,*) '  NSPEC             : ',MGLL_QMU_V%nspec
    write(IMAIN,*) '  NGLLX/NGLLY/NGLLZ : ',NGLLX,NGLLY,NGLLZ
    call flush_IMAIN()
  endif

  ! safety check
  if (MGLL_QMU_V%nspec /= MGLL_V%nspec) &
    stop 'Invalid nspec for QMu GLL model, size must be same to model GLL'

  ! allocates arrays
  ! differs for isotropic model or transverse isotropic models
  allocate( MGLL_QMU_V%qmu_new(NGLLX,NGLLY,NGLLZ,MGLL_QMU_V%nspec),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating qmu_new array')
  MGLL_QMU_V%qmu_new(:,:,:,:) = 0.0_CUSTOM_REAL

  ! reads in model files for each process
  if (ADIOS_FOR_MODELS) then
     call read_gll_qmu_model_adios(myrank)
  else
     call read_gll_qmu_model(myrank)
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'model Qmu GLL:'
    call flush_IMAIN()
  endif
  call print_gll_min_max_all(MGLL_QMU_V%nspec,MGLL_QMU_V%qmu_new,"qmu new")

  call synchronize_all()

  end subroutine model_attenuation_gll_broadcast

!
!-------------------------------------------------------------------------------------------------
!


  subroutine read_gll_qmu_model(rank)

  use constants
  use model_gll_qmu_par

  implicit none

  integer,intent(in) :: rank

  ! local parameters
  integer :: ier
  character(len=MAX_STRING_LEN) :: filename
  character(len=MAX_STRING_LEN) :: prname

  if (myrank == 0) then
    write(IMAIN,*) 'reading in model from: ',trim(PATHNAME_GLL_modeldir)
    if (rank /= myrank) write(IMAIN,*) '  mesh slice for rank: ',rank
    call flush_IMAIN()
  endif

  ! for example: DATA/GLL/proc000000_reg1_qmu.bin
  !
  ! root name
  write(prname,'(a,i6.6,a)') PATHNAME_GLL_modeldir(1:len_trim(PATHNAME_GLL_modeldir))//'proc',rank,'_reg1_'

  ! qmu mesh
  filename = prname(1:len_trim(prname))//'qmu.bin'
  open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(filename)
    call exit_MPI(rank,'Error model Qmu GLL')
  endif
  read(IIN) MGLL_QMU_V%qmu_new(:,:,:,1:MGLL_QMU_V%nspec)
  close(IIN)

  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  reading done'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine read_gll_qmu_model


!
!-------------------------------------------------------------------------------------------------
!


  subroutine model_attenuation_gll(ispec, i, j, k, Qmu)

  use model_gll_qmu_par

  implicit none

  integer,intent(in) :: ispec,i,j,k
  double precision,intent(inout) :: Qmu

  Qmu = dble( MGLL_QMU_V%qmu_new(i,j,k,ispec) )

  end subroutine model_attenuation_gll
