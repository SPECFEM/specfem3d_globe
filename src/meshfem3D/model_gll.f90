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

!--------------------------------------------------------------------------------------------------
! GLL
!
! based on modified GLL mesh output from mesher
!
! used for iterative inversion procedures
!--------------------------------------------------------------------------------------------------

  subroutine model_gll_broadcast(myrank,MGLL_V,NSPEC)

! standard routine to setup model

  use constants
  use meshfem3D_models_par,only: TRANSVERSE_ISOTROPY
  use meshfem3D_par, only: ADIOS_FOR_MODELS

  implicit none

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
  double precision :: scaleval
  real(kind=CUSTOM_REAL) :: minvalue,maxvalue,min_all,max_all
  integer :: ier

  ! allocates arrays
  ! differs for isotropic model or transverse isotropic models
  if (.not. TRANSVERSE_ISOTROPY) then
    ! isotropic model
    allocate( MGLL_V%vp_new(NGLLX,NGLLY,NGLLZ,NSPEC(IREGION_CRUST_MANTLE)), &
              MGLL_V%vs_new(NGLLX,NGLLY,NGLLZ,NSPEC(IREGION_CRUST_MANTLE)), stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating vp_new,.. arrays')
  else
    ! transverse isotropic model
    allocate( MGLL_V%vpv_new(NGLLX,NGLLY,NGLLZ,NSPEC(IREGION_CRUST_MANTLE)), &
              MGLL_V%vph_new(NGLLX,NGLLY,NGLLZ,NSPEC(IREGION_CRUST_MANTLE)), &
              MGLL_V%vsv_new(NGLLX,NGLLY,NGLLZ,NSPEC(IREGION_CRUST_MANTLE)), &
              MGLL_V%vsh_new(NGLLX,NGLLY,NGLLZ,NSPEC(IREGION_CRUST_MANTLE)), &
              MGLL_V%eta_new(NGLLX,NGLLY,NGLLZ,NSPEC(IREGION_CRUST_MANTLE)), stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating vpv_new,.. arrays')

  endif
  allocate( MGLL_V%rho_new(NGLLX,NGLLY,NGLLZ,NSPEC(IREGION_CRUST_MANTLE)), stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating rho_new,.. arrays')

  ! reads in model files for each process
  if (ADIOS_FOR_MODELS) then
    call read_gll_model_adios(myrank,MGLL_V,NSPEC)
  else
    call read_gll_model(myrank,MGLL_V,NSPEC)
  endif

  ! checks velocity range
  if (.not. TRANSVERSE_ISOTROPY) then

    ! isotropic model
    if (myrank == 0) then
      write(IMAIN,*)'model GLL: isotropic'
      call flush_IMAIN()
    endif

    ! Vs
    maxvalue = maxval( MGLL_V%vs_new )
    minvalue = minval( MGLL_V%vs_new )
    call max_all_cr(maxvalue, max_all)
    call min_all_cr(minvalue, min_all)
    if (myrank == 0) then
      write(IMAIN,*) '  vs new min/max: ',min_all,max_all
    endif
    ! Vp
    maxvalue = maxval( MGLL_V%vp_new )
    minvalue = minval( MGLL_V%vp_new )
    call max_all_cr(maxvalue, max_all)
    call min_all_cr(minvalue, min_all)
    if (myrank == 0) then
      write(IMAIN,*) '  vp new min/max: ',min_all,max_all
    endif
    ! density
    maxvalue = maxval( MGLL_V%rho_new )
    minvalue = minval( MGLL_V%rho_new )
    call max_all_cr(maxvalue, max_all)
    call min_all_cr(minvalue, min_all)
    if (myrank == 0) then
      write(IMAIN,*) '  rho new min/max: ',min_all,max_all
      write(IMAIN,*)
      call flush_IMAIN()
    endif

  else

    ! transverse isotropic model
    if (myrank == 0) then
      write(IMAIN,*)'model GLL: transverse isotropic'
    endif

    ! Vsv
    maxvalue = maxval( MGLL_V%vsv_new )
    minvalue = minval( MGLL_V%vsv_new )
    call max_all_cr(maxvalue, max_all)
    call min_all_cr(minvalue, min_all)
    if (myrank == 0) then
      write(IMAIN,*) '  vsv new min/max: ',min_all,max_all
    endif
    ! Vsh
    maxvalue = maxval( MGLL_V%vsh_new )
    minvalue = minval( MGLL_V%vsh_new )
    call max_all_cr(maxvalue, max_all)
    call min_all_cr(minvalue, min_all)
    if (myrank == 0) then
      write(IMAIN,*) '  vsh new min/max: ',min_all,max_all
    endif
    ! Vpv
    maxvalue = maxval( MGLL_V%vpv_new )
    minvalue = minval( MGLL_V%vpv_new )
    call max_all_cr(maxvalue, max_all)
    call min_all_cr(minvalue, min_all)
    if (myrank == 0) then
      write(IMAIN,*) '  vpv new min/max: ',min_all,max_all
    endif
    ! Vph
    maxvalue = maxval( MGLL_V%vph_new )
    minvalue = minval( MGLL_V%vph_new )
    call max_all_cr(maxvalue, max_all)
    call min_all_cr(minvalue, min_all)
    if (myrank == 0) then
      write(IMAIN,*) '  vph new min/max: ',min_all,max_all
    endif
    ! density
    maxvalue = maxval( MGLL_V%rho_new )
    minvalue = minval( MGLL_V%rho_new )
    call max_all_cr(maxvalue, max_all)
    call min_all_cr(minvalue, min_all)
    if (myrank == 0) then
      write(IMAIN,*) '  rho new min/max: ',min_all,max_all
    endif
    ! eta
    maxvalue = maxval( MGLL_V%eta_new )
    minvalue = minval( MGLL_V%eta_new )
    call max_all_cr(maxvalue, max_all)
    call min_all_cr(minvalue, min_all)
    if (myrank == 0) then
      write(IMAIN,*) '  eta new min/max: ',min_all,max_all
      write(IMAIN,*)
      call flush_IMAIN()
    endif

  endif

  ! non-dimensionalizes model values
  ! (SPECFEM3D_GLOBE uses non-dimensionalized values in subsequent computations)
  ! scaling values
  ! (model velocities must be given as km/s)
  scaleval = dsqrt(PI*GRAV*RHOAV)
  MGLL_V%scale_velocity = 1000.0d0/(R_EARTH*scaleval)
  MGLL_V%scale_density =  1000.0d0/RHOAV
  if (.not. TRANSVERSE_ISOTROPY) then
      ! non-dimensionalize isotropic values
      MGLL_V%vp_new = MGLL_V%vp_new * MGLL_V%scale_velocity
      MGLL_V%vs_new = MGLL_V%vs_new * MGLL_V%scale_velocity
      MGLL_V%rho_new = MGLL_V%rho_new * MGLL_V%scale_density
  else
      ! non-dimensionalize
      ! transverse isotropic model
      MGLL_V%vpv_new = MGLL_V%vpv_new * MGLL_V%scale_velocity
      MGLL_V%vph_new = MGLL_V%vph_new * MGLL_V%scale_velocity
      MGLL_V%vsv_new = MGLL_V%vsv_new * MGLL_V%scale_velocity
      MGLL_V%vsh_new = MGLL_V%vsh_new * MGLL_V%scale_velocity
      MGLL_V%rho_new = MGLL_V%rho_new * MGLL_V%scale_density
      ! eta is already non-dimensional
  endif

  end subroutine model_gll_broadcast

!
!-------------------------------------------------------------------------------------------------
!


  subroutine read_gll_model(myrank,MGLL_V,NSPEC)

  use constants
  use meshfem3D_models_par,only: TRANSVERSE_ISOTROPY

  implicit none

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
  integer :: ier
  character(len=MAX_STRING_LEN) :: prname

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*)'reading in model from ',trim(PATHNAME_GLL_modeldir)
    call flush_IMAIN()
  endif

  ! only crust and mantle
  write(prname,'(a,i6.6,a)') PATHNAME_GLL_modeldir(1:len_trim(PATHNAME_GLL_modeldir))//'proc',myrank,'_reg1_'

  ! reads in model for each partition
  if (.not. TRANSVERSE_ISOTROPY) then
    ! isotropic model
    ! vp mesh
    open(unit=IIN,file=prname(1:len_trim(prname))//'vp.bin',&
          status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      write(IMAIN,*) 'Error opening: ',prname(1:len_trim(prname))//'vp.bin'
      call exit_MPI(myrank,'Error model GLL')
    endif
    read(IIN) MGLL_V%vp_new(:,:,:,1:nspec(IREGION_CRUST_MANTLE))
    close(IIN)

    ! vs mesh
    open(unit=IIN,file=prname(1:len_trim(prname))//'vs.bin', &
         status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print*,'Error opening: ',prname(1:len_trim(prname))//'vs.bin'
      call exit_MPI(myrank,'Error model GLL')
    endif
    read(IIN) MGLL_V%vs_new(:,:,:,1:nspec(IREGION_CRUST_MANTLE))
    close(IIN)

  else

    ! transverse isotropic model
    ! vp mesh
    open(unit=IIN,file=prname(1:len_trim(prname))//'vpv.bin',&
          status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      write(IMAIN,*) 'Error opening: ',prname(1:len_trim(prname))//'vpv.bin'
      call exit_MPI(myrank,'Error model GLL')
    endif
    read(IIN) MGLL_V%vpv_new(:,:,:,1:nspec(IREGION_CRUST_MANTLE))
    close(IIN)

    open(unit=IIN,file=prname(1:len_trim(prname))//'vph.bin',&
          status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      write(IMAIN,*) 'Error opening: ',prname(1:len_trim(prname))//'vph.bin'
      call exit_MPI(myrank,'Error model GLL')
    endif
    read(IIN) MGLL_V%vph_new(:,:,:,1:nspec(IREGION_CRUST_MANTLE))
    close(IIN)

    ! vs mesh
    open(unit=IIN,file=prname(1:len_trim(prname))//'vsv.bin', &
         status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print*,'Error opening: ',prname(1:len_trim(prname))//'vsv.bin'
      call exit_MPI(myrank,'Error model GLL')
    endif
    read(IIN) MGLL_V%vsv_new(:,:,:,1:nspec(IREGION_CRUST_MANTLE))
    close(IIN)

    open(unit=IIN,file=prname(1:len_trim(prname))//'vsh.bin', &
         status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print*,'Error opening: ',prname(1:len_trim(prname))//'vsh.bin'
      call exit_MPI(myrank,'Error model GLL')
    endif
    read(IIN) MGLL_V%vsh_new(:,:,:,1:nspec(IREGION_CRUST_MANTLE))
    close(IIN)

    ! eta mesh
    open(unit=IIN,file=prname(1:len_trim(prname))//'eta.bin', &
         status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print*,'Error opening: ',prname(1:len_trim(prname))//'eta.bin'
      call exit_MPI(myrank,'Error model GLL')
    endif
    read(IIN) MGLL_V%eta_new(:,:,:,1:nspec(IREGION_CRUST_MANTLE))
    close(IIN)

  endif

  ! rho mesh
  open(unit=IIN,file=prname(1:len_trim(prname))//'rho.bin', &
       status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print*,'Error opening: ',prname(1:len_trim(prname))//'rho.bin'
    call exit_MPI(myrank,'Error model GLL')
  endif
  read(IIN) MGLL_V%rho_new(:,:,:,1:nspec(IREGION_CRUST_MANTLE))
  close(IIN)

  end subroutine read_gll_model
