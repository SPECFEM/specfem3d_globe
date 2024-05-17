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
! GLL
!
! based on modified GLL mesh output from mesher
!
! used for iterative inversion procedures
!--------------------------------------------------------------------------------------------------

  module model_gll_par

  use constants, only: CUSTOM_REAL

  ! GLL model_variables
  type model_gll_variables
    sequence
    ! tomographic iteration model on GLL points
    real(kind=CUSTOM_REAL) :: scale_velocity,scale_density,scale_GPa

    ! isotropic model
    real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: vs_new,vp_new,rho_new

    ! transverse isotropic model
    real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: vsv_new,vpv_new, &
                                                             vsh_new,vph_new,eta_new
    ! azimuthal anisotropy
    real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: Gc_prime_new,Gs_prime_new,mu0_new

    ! number of elements (crust/mantle elements)
    integer :: nspec
    integer :: dummy_pad ! padding 4 bytes to align the structure
  end type model_gll_variables

  ! crust/mantle model
  type (model_gll_variables) :: MGLL_V

  ! inner core model
  type (model_gll_variables) :: MGLL_V_IC

  ! model GLL type: 1 == iso, 2 == tiso, 3 == azi
  integer :: MGLL_TYPE

  end module model_gll_par

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_gll_broadcast()

! standard routine to setup model

  use constants
  use shared_parameters, only: NSPEC_REGIONS,ADIOS_FOR_MODELS,NPROCTOT,NCHUNKS, &
    MODEL,MODEL_GLL_TYPE,R_PLANET,RHOAV

  use model_gll_par

  implicit none

  ! local parameters
  double precision :: scaleval
  integer :: ier

  ! sets type (iso,tiso,azi)
  MGLL_TYPE = MODEL_GLL_TYPE

  ! sets number of elements (crust/mantle region)
  MGLL_V%nspec = NSPEC_REGIONS(IREGION_CRUST_MANTLE)
  MGLL_V_IC%nspec = 0

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

  ! saftey check
  if (MGLL_TYPE < 1 .or. MGLL_TYPE > 3) &
    stop 'Invalid MODEL_GLL_TYPE, please use 1(iso), 2(tiso) or 3(azi) in get_model_parameters.F90 setup'

  ! allocates arrays
  ! differs for isotropic model or transverse isotropic models
  select case(MGLL_TYPE)
  case (1)
    ! isotropic model
    allocate( MGLL_V%vp_new(NGLLX,NGLLY,NGLLZ,MGLL_V%nspec), &
              MGLL_V%vs_new(NGLLX,NGLLY,NGLLZ,MGLL_V%nspec), stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating vp_new,.. arrays')
    MGLL_V%vp_new(:,:,:,:) = 0.0_CUSTOM_REAL
    MGLL_V%vs_new(:,:,:,:) = 0.0_CUSTOM_REAL
  case (2)
    ! transverse isotropic model
    allocate( MGLL_V%vpv_new(NGLLX,NGLLY,NGLLZ,MGLL_V%nspec), &
              MGLL_V%vph_new(NGLLX,NGLLY,NGLLZ,MGLL_V%nspec), &
              MGLL_V%vsv_new(NGLLX,NGLLY,NGLLZ,MGLL_V%nspec), &
              MGLL_V%vsh_new(NGLLX,NGLLY,NGLLZ,MGLL_V%nspec), &
              MGLL_V%eta_new(NGLLX,NGLLY,NGLLZ,MGLL_V%nspec), stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating vpv_new,.. arrays')
    MGLL_V%vpv_new(:,:,:,:) = 0.0_CUSTOM_REAL
    MGLL_V%vph_new(:,:,:,:) = 0.0_CUSTOM_REAL
    MGLL_V%vsv_new(:,:,:,:) = 0.0_CUSTOM_REAL
    MGLL_V%vsh_new(:,:,:,:) = 0.0_CUSTOM_REAL
    MGLL_V%eta_new(:,:,:,:) = 0.0_CUSTOM_REAL
  case (3)
    ! azimuthally anisotropic model
    allocate( MGLL_V%vpv_new(NGLLX,NGLLY,NGLLZ,MGLL_V%nspec), &
              MGLL_V%vph_new(NGLLX,NGLLY,NGLLZ,MGLL_V%nspec), &
              MGLL_V%vsv_new(NGLLX,NGLLY,NGLLZ,MGLL_V%nspec), &
              MGLL_V%vsh_new(NGLLX,NGLLY,NGLLZ,MGLL_V%nspec), &
              MGLL_V%eta_new(NGLLX,NGLLY,NGLLZ,MGLL_V%nspec), stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating vpv_new,.. arrays')
    MGLL_V%vpv_new(:,:,:,:) = 0.0_CUSTOM_REAL
    MGLL_V%vph_new(:,:,:,:) = 0.0_CUSTOM_REAL
    MGLL_V%vsv_new(:,:,:,:) = 0.0_CUSTOM_REAL
    MGLL_V%vsh_new(:,:,:,:) = 0.0_CUSTOM_REAL
    MGLL_V%eta_new(:,:,:,:) = 0.0_CUSTOM_REAL
    allocate( MGLL_V%Gc_prime_new(NGLLX,NGLLY,NGLLZ,MGLL_V%nspec), &
              MGLL_V%Gs_prime_new(NGLLX,NGLLY,NGLLZ,MGLL_V%nspec), &
              MGLL_V%mu0_new(NGLLX,NGLLY,NGLLZ,MGLL_V%nspec),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating vpv_new,.. arrays')
    MGLL_V%Gc_prime_new(:,:,:,:) = 0.0_CUSTOM_REAL
    MGLL_V%Gs_prime_new(:,:,:,:) = 0.0_CUSTOM_REAL
    MGLL_V%mu0_new(:,:,:,:) = 0.0_CUSTOM_REAL
  case default
    stop 'Invalid MGLL_TYPE, type not implemented yet'
  end select

  allocate( MGLL_V%rho_new(NGLLX,NGLLY,NGLLZ,MGLL_V%nspec), stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating rho_new,.. arrays')
  MGLL_V%rho_new(:,:,:,:) = 0.0_CUSTOM_REAL

  ! reads in model files for each process
  if (ADIOS_FOR_MODELS) then
    call read_gll_model_adios(myrank)
  else
    call read_gll_model(myrank)
  endif

  ! outputs velocity range
  select case(MGLL_TYPE)
  case (1)
    ! isotropic model
    if (myrank == 0) then
      write(IMAIN,*) 'model GLL: isotropic'
      call flush_IMAIN()
    endif
    ! Vp
    call print_gll_min_max_all(MGLL_V%nspec,MGLL_V%vp_new,"vp new")
    ! Vs
    call print_gll_min_max_all(MGLL_V%nspec,MGLL_V%vs_new,"vs new")

  case (2)
    ! transverse
    if (myrank == 0) then
      ! transverse isotropic model
      write(IMAIN,*) 'model GLL: transverse isotropic'
      call flush_IMAIN()
    endif
    ! Vpv
    call print_gll_min_max_all(MGLL_V%nspec,MGLL_V%vpv_new,"vpv new")
    ! Vph
    call print_gll_min_max_all(MGLL_V%nspec,MGLL_V%vph_new,"vph new")
    ! Vsv
    call print_gll_min_max_all(MGLL_V%nspec,MGLL_V%vsv_new,"vsv new")
    ! Vsh
    call print_gll_min_max_all(MGLL_V%nspec,MGLL_V%vsh_new,"vsh new")
    ! eta
    call print_gll_min_max_all(MGLL_V%nspec,MGLL_V%eta_new,"eta new")

  case (3)
    ! azimuthal model
    if (myrank == 0) then
      ! azimuthal anisotropic model
      write(IMAIN,*) 'model GLL: azimuthal anisotropic'
      call flush_IMAIN()
    endif
    ! Vpv
    call print_gll_min_max_all(MGLL_V%nspec,MGLL_V%vpv_new,"vpv new")
    ! Vph
    call print_gll_min_max_all(MGLL_V%nspec,MGLL_V%vph_new,"vph new")
    ! Vsv
    call print_gll_min_max_all(MGLL_V%nspec,MGLL_V%vsv_new,"vsv new")
    ! Vsh
    call print_gll_min_max_all(MGLL_V%nspec,MGLL_V%vsh_new,"vsh new")
    ! eta
    call print_gll_min_max_all(MGLL_V%nspec,MGLL_V%eta_new,"eta new")
    ! Gc_prime
    call print_gll_min_max_all(MGLL_V%nspec,MGLL_V%Gc_prime_new,"Gc_prime new")
    ! Gs_prime
    call print_gll_min_max_all(MGLL_V%nspec,MGLL_V%Gs_prime_new,"Gs_prime new")
    ! mu0
    call print_gll_min_max_all(MGLL_V%nspec,MGLL_V%mu0_new,"mu0 new")

  case default
    stop 'New MGLL_TYPE, velocity range check not implemented yet'
  end select

  ! density
  call print_gll_min_max_all(MGLL_V%nspec,MGLL_V%rho_new,"rho new")
  if (myrank == 0) then
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! non-dimensionalizes model values
  ! (SPECFEM3D_GLOBE uses non-dimensionalized values in subsequent computations)
  ! scaling values
  ! (model velocities must be given as km/s)
  scaleval = dsqrt(PI*GRAV*RHOAV)
  MGLL_V%scale_velocity = real(1000.0d0/(R_PLANET*scaleval),kind=CUSTOM_REAL)
  MGLL_V%scale_density  = real(1000.0d0/RHOAV,kind=CUSTOM_REAL)
  ! non-dimensionalize the elastic coefficients using the scale of GPa--[g/cm^3][(km/s)^2]
  ! equal to scale_density * scale_velocity**2
  MGLL_V%scale_GPa      = real(1.d0/( (RHOAV/1000.d0)*((R_PLANET*scaleval/1000.d0)**2) ),kind=CUSTOM_REAL)

  select case(MGLL_TYPE)
  case (1)
    ! non-dimensionalize isotropic values
    MGLL_V%vp_new = MGLL_V%vp_new * MGLL_V%scale_velocity
    MGLL_V%vs_new = MGLL_V%vs_new * MGLL_V%scale_velocity
    MGLL_V%rho_new = MGLL_V%rho_new * MGLL_V%scale_density
  case (2)
    ! non-dimensionalize
    ! transverse isotropic model
    MGLL_V%vpv_new = MGLL_V%vpv_new * MGLL_V%scale_velocity
    MGLL_V%vph_new = MGLL_V%vph_new * MGLL_V%scale_velocity
    MGLL_V%vsv_new = MGLL_V%vsv_new * MGLL_V%scale_velocity
    MGLL_V%vsh_new = MGLL_V%vsh_new * MGLL_V%scale_velocity
    MGLL_V%rho_new = MGLL_V%rho_new * MGLL_V%scale_density
    ! eta is already non-dimensional
  case (3)
    ! azimuthal model
    MGLL_V%vpv_new = MGLL_V%vpv_new * MGLL_V%scale_velocity
    MGLL_V%vph_new = MGLL_V%vph_new * MGLL_V%scale_velocity
    MGLL_V%vsv_new = MGLL_V%vsv_new * MGLL_V%scale_velocity
    MGLL_V%vsh_new = MGLL_V%vsh_new * MGLL_V%scale_velocity
    MGLL_V%rho_new = MGLL_V%rho_new * MGLL_V%scale_density
    ! Gc_prime,Gs_prime are already non-dimensional
    MGLL_V%mu0_new = MGLL_V%mu0_new * MGLL_V%scale_GPa ! moduli in GPa

    ! eta is already non-dimensional
  end select

  ! (optional) inner core
  if (MGLL_V_IC%nspec > 0) then
    ! only isotropic
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  with inner core (vp,vs,rho):'
      call flush_IMAIN()
    endif
    call print_gll_min_max_all(MGLL_V_IC%nspec,MGLL_V_IC%vp_new,"vp new")
    call print_gll_min_max_all(MGLL_V_IC%nspec,MGLL_V_IC%vs_new,"vs new")
    call print_gll_min_max_all(MGLL_V_IC%nspec,MGLL_V_IC%rho_new,"rho new")
    if (myrank == 0) then
      write(IMAIN,*)
      call flush_IMAIN()
    endif
    ! non-dimensionalizes
    MGLL_V_IC%scale_velocity = real(1000.0d0/(R_PLANET*scaleval),kind=CUSTOM_REAL)
    MGLL_V_IC%scale_density  = real(1000.0d0/RHOAV,kind=CUSTOM_REAL)

    MGLL_V_IC%vp_new = MGLL_V_IC%vp_new * MGLL_V_IC%scale_velocity
    MGLL_V_IC%vs_new = MGLL_V_IC%vs_new * MGLL_V_IC%scale_velocity
    MGLL_V_IC%rho_new = MGLL_V_IC%rho_new * MGLL_V_IC%scale_density
  endif

  call synchronize_all()

  end subroutine model_gll_broadcast

!
!-------------------------------------------------------------------------------------------------
!


  subroutine read_gll_model(rank)

  use constants, only: MAX_STRING_LEN,IMAIN,NGLLX,NGLLY,NGLLZ,PATHNAME_GLL_modeldir,myrank

  use shared_parameters, only: NCHUNKS,NPROCTOT,NPROC_XI,NPROC_ETA,NEX_XI,NEX_ETA

  use model_gll_par

  implicit none

  integer,intent(in) :: rank

  ! local parameters
  integer :: num_ranks,nproc_eta_gll,nproc_xi_gll,nspec_gll,nex_xi_gll,nex_eta_gll
  integer(kind=8) :: filesize
  logical :: exists
  character(len=MAX_STRING_LEN) :: filename

  if (myrank == 0) then
    write(IMAIN,*) 'reading in model from: ',trim(PATHNAME_GLL_modeldir)
    if (rank /= myrank) write(IMAIN,*) '  mesh slice for rank: ',rank
    call flush_IMAIN()
  endif

  ! gets number of processes from the mesh used for this GLL model (1 file per process)
  if (myrank == 0) then
    ! counts number of files and estimates setup
    num_ranks = 0
    filesize = 0
    exists = .true.
    do while (exists)
      ! checks with density rho filenames
      write(filename,'(a,i6.6,a)') PATHNAME_GLL_modeldir(1:len_trim(PATHNAME_GLL_modeldir))//'proc',num_ranks,'_reg1_rho.bin'
      inquire(file=trim(filename),exist=exists)
      if (exists) then
        ! gets file size in bytes
        if (num_ranks == 0) then
          inquire(file=trim(filename),size=filesize)
        endif
        num_ranks = num_ranks + 1
      endif
    enddo
    ! checks
    if (num_ranks == 0) then
      print *,'Error invalid number of rho mesh-files found for GLL model'
      stop 'Error invalid number of rho mesh-files found for GLL model'
    endif
    ! assumes same number of chunks and nproc_eta = nproc_xi
    nproc_eta_gll = int (sqrt ( dble(num_ranks) / NCHUNKS ))
    nproc_xi_gll = int (sqrt ( dble(num_ranks) / NCHUNKS ))

    ! estimates nspec (filesize might have 2byte overhead)
    nspec_gll = int (filesize / dble(NGLLX*NGLLY*NGLLZ) / dble(CUSTOM_REAL) )
    nex_xi_gll = int (sqrt(dble(nspec_gll) / dble(MGLL_V%nspec)) * NEX_XI )
    nex_eta_gll = int (sqrt(dble(nspec_gll) / dble(MGLL_V%nspec)) * NEX_ETA )

    ! user info
    write(IMAIN,*) '  file setting found:'
    write(IMAIN,*) '                    filesize : ',filesize
    write(IMAIN,*) '    total number of processes: ',num_ranks
    write(IMAIN,*) '                     NPROC_XI: ',nproc_xi_gll
    write(IMAIN,*) '                    NPROC_ETA: ',nproc_eta_gll
    write(IMAIN,*) '              estimated nspec: ',nspec_gll
    write(IMAIN,*) '     estimated nex_xi/nex_eta: ',nex_xi_gll,nex_eta_gll
    write(IMAIN,*)
    write(IMAIN,*) '  current setting:'
    write(IMAIN,*) '    total number of processes: ',NPROCTOT
    write(IMAIN,*) '                     NPROC_XI: ',NPROC_XI
    write(IMAIN,*) '                    NPROC_ETA: ',NPROC_ETA
    write(IMAIN,*) '                        nspec: ',MGLL_V%nspec
    write(IMAIN,*) '               nex_xi/nex_eta: ',NEX_XI,NEX_ETA
    write(IMAIN,*)
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! only crust and mantle has GLL model files for now
  call read_gll_modelfiles(rank)

  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  reading done'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine read_gll_model


!
!-------------------------------------------------------------------------------------------------
!


  subroutine read_gll_modelfiles(rank)

  use constants, only: MAX_STRING_LEN,IMAIN,IIN,PATHNAME_GLL_modeldir, &
    NGLLX,NGLLY,NGLLZ,IREGION_INNER_CORE

  use shared_parameters, only: NSPEC_REGIONS

  use model_gll_par

  implicit none

  integer,intent(in) :: rank

  ! local parameters
  integer :: ier
  character(len=MAX_STRING_LEN) :: prname,filename
  logical :: has_innercore,has_innercore_all

  ! crust/mantle model
  ! for example: DATA/GLL/proc000000_reg1_rho.bin
  !
  ! root name
  write(prname,'(a,i6.6,a)') PATHNAME_GLL_modeldir(1:len_trim(PATHNAME_GLL_modeldir))//'proc',rank,'_reg1_'

  ! reads in model for each partition
  select case (MGLL_TYPE)
  case (1)
    ! isotropic model
    if (rank == 0) then
      write(IMAIN,*) '  reads isotropic model values: rho,vp,vs'
      call flush_IMAIN()
    endif

    ! vp mesh
    filename = prname(1:len_trim(prname))//'vp.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      write(IMAIN,*) 'Error opening: ',trim(filename)
      call exit_MPI(rank,'Error model GLL')
    endif
    read(IIN) MGLL_V%vp_new(:,:,:,1:MGLL_V%nspec)
    close(IIN)

    ! vs mesh
    filename = prname(1:len_trim(prname))//'vs.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',trim(filename)
      call exit_MPI(rank,'Error model GLL')
    endif
    read(IIN) MGLL_V%vs_new(:,:,:,1:MGLL_V%nspec)
    close(IIN)

  case (2)
    ! transverse isotropic model
    if (rank == 0) then
      write(IMAIN,*) '  reads transversely isotropic model values: rho,vpv,vph,vsv,vsh,eta'
      call flush_IMAIN()
    endif

    ! vpv/vph mesh
    filename = prname(1:len_trim(prname))//'vpv.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      write(IMAIN,*) 'Error opening: ',trim(filename)
      call exit_MPI(rank,'Error model GLL')
    endif
    read(IIN) MGLL_V%vpv_new(:,:,:,1:MGLL_V%nspec)
    close(IIN)

    filename = prname(1:len_trim(prname))//'vph.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      write(IMAIN,*) 'Error opening: ',trim(filename)
      call exit_MPI(rank,'Error model GLL')
    endif
    read(IIN) MGLL_V%vph_new(:,:,:,1:MGLL_V%nspec)
    close(IIN)

    ! vsv/vsh mesh
    filename = prname(1:len_trim(prname))//'vsv.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',trim(filename)
      call exit_MPI(rank,'Error model GLL')
    endif
    read(IIN) MGLL_V%vsv_new(:,:,:,1:MGLL_V%nspec)
    close(IIN)

    filename = prname(1:len_trim(prname))//'vsh.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',trim(filename)
      call exit_MPI(rank,'Error model GLL')
    endif
    read(IIN) MGLL_V%vsh_new(:,:,:,1:MGLL_V%nspec)
    close(IIN)

    ! eta mesh
    filename = prname(1:len_trim(prname))//'eta.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',trim(filename)
      call exit_MPI(rank,'Error model GLL')
    endif
    read(IIN) MGLL_V%eta_new(:,:,:,1:MGLL_V%nspec)
    close(IIN)

  case (3)
    ! azimuthal model
    if (rank == 0) then
      write(IMAIN,*) '  reads azimuthal anisotropic model values: rho,vpv,vph,vsv,vsh,eta,Gc_prime,Gs_prime,mu0'
      call flush_IMAIN()
    endif

    ! vpv/vph mesh
    filename = prname(1:len_trim(prname))//'vpv.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      write(IMAIN,*) 'Error opening: ',trim(filename)
      call exit_MPI(rank,'Error model GLL')
    endif
    read(IIN) MGLL_V%vpv_new(:,:,:,1:MGLL_V%nspec)
    close(IIN)

    filename = prname(1:len_trim(prname))//'vph.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      write(IMAIN,*) 'Error opening: ',trim(filename)
      call exit_MPI(rank,'Error model GLL')
    endif
    read(IIN) MGLL_V%vph_new(:,:,:,1:MGLL_V%nspec)
    close(IIN)

    ! vsv/vsh mesh
    filename = prname(1:len_trim(prname))//'vsv.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',trim(filename)
      call exit_MPI(rank,'Error model GLL')
    endif
    read(IIN) MGLL_V%vsv_new(:,:,:,1:MGLL_V%nspec)
    close(IIN)

    filename = prname(1:len_trim(prname))//'vsh.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',trim(filename)
      call exit_MPI(rank,'Error model GLL')
    endif
    read(IIN) MGLL_V%vsh_new(:,:,:,1:MGLL_V%nspec)
    close(IIN)

    ! eta mesh
    filename = prname(1:len_trim(prname))//'eta.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',trim(filename)
      call exit_MPI(rank,'Error model GLL')
    endif
    read(IIN) MGLL_V%eta_new(:,:,:,1:MGLL_V%nspec)
    close(IIN)

    ! Gc_prime
    filename = prname(1:len_trim(prname))//'Gc_prime.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',trim(filename)
      call exit_MPI(rank,'Error model GLL')
    endif
    read(IIN) MGLL_V%Gc_prime_new(:,:,:,1:MGLL_V%nspec)
    close(IIN)

    ! Gs_prime
    filename = prname(1:len_trim(prname))//'Gs_prime.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',trim(filename)
      call exit_MPI(rank,'Error model GLL')
    endif
    read(IIN) MGLL_V%Gs_prime_new(:,:,:,1:MGLL_V%nspec)
    close(IIN)

    ! mu0
    if (.true.) then  ! by default try to read in mu0
      filename = prname(1:len_trim(prname))//'mu0.bin'
      open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
      if (ier /= 0) then
        print *,'Error opening: ',trim(filename)
        call exit_MPI(rank,'Error model GLL')
      endif
      read(IIN) MGLL_V%mu0_new(:,:,:,1:MGLL_V%nspec)
      close(IIN)
    else
      MGLL_V%mu0_new(:,:,:,1:MGLL_V%nspec) = 1.0
    endif

  case default
    stop 'Invalid MGLL_TYPE, mesh file reading not implemented yet'
  end select

  ! rho mesh
  filename = prname(1:len_trim(prname))//'rho.bin'
  open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(filename)
    call exit_MPI(rank,'Error model GLL')
  endif
  read(IIN) MGLL_V%rho_new(:,:,:,1:MGLL_V%nspec)
  close(IIN)

  ! inner core model
  ! for example: DATA/GLL/proc000000_reg3_rho.bin
  !
  ! root name
  write(prname,'(a,i6.6,a)') PATHNAME_GLL_modeldir(1:len_trim(PATHNAME_GLL_modeldir))//'proc',rank,'_reg3_'

  ! checks if rho file exists for inner core
  filename = prname(1:len_trim(prname))//'rho.bin'
  inquire(file=trim(filename),exist=has_innercore)

  ! if file found, all processes will try to read in inner core files
  call any_all_l(has_innercore,has_innercore_all)

  ! reads in files
  if (has_innercore_all) then
    ! only isotropic model file supported so far for inner core
    MGLL_V_IC%nspec = NSPEC_REGIONS(IREGION_INNER_CORE)

    ! allocates arrays
    allocate( MGLL_V_IC%vp_new(NGLLX,NGLLY,NGLLZ,MGLL_V_IC%nspec), &
              MGLL_V_IC%vs_new(NGLLX,NGLLY,NGLLZ,MGLL_V_IC%nspec), stat=ier)
    if (ier /= 0 ) call exit_MPI(rank,'Error allocating inner core vp_new,.. arrays')
    MGLL_V_IC%vp_new(:,:,:,:) = 0.0_CUSTOM_REAL
    MGLL_V_IC%vs_new(:,:,:,:) = 0.0_CUSTOM_REAL

    allocate( MGLL_V_IC%rho_new(NGLLX,NGLLY,NGLLZ,MGLL_V_IC%nspec), stat=ier)
    if (ier /= 0 ) call exit_MPI(rank,'Error allocating inner core rho_new,.. arrays')
    MGLL_V_IC%rho_new(:,:,:,:) = 0.0_CUSTOM_REAL

    ! rho mesh
    filename = prname(1:len_trim(prname))//'rho.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',trim(filename)
      call exit_MPI(rank,'Error model GLL')
    endif
    read(IIN) MGLL_V_IC%rho_new(:,:,:,1:MGLL_V_IC%nspec)
    close(IIN)
    ! vp mesh
    filename = prname(1:len_trim(prname))//'vp.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      write(IMAIN,*) 'Error opening: ',trim(filename)
      call exit_MPI(rank,'Error model GLL')
    endif
    read(IIN) MGLL_V_IC%vp_new(:,:,:,1:MGLL_V_IC%nspec)
    close(IIN)
    ! vs mesh
    filename = prname(1:len_trim(prname))//'vs.bin'
    open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',trim(filename)
      call exit_MPI(rank,'Error model GLL')
    endif
    read(IIN) MGLL_V_IC%vs_new(:,:,:,1:MGLL_V_IC%nspec)
    close(IIN)
  endif

  end subroutine read_gll_modelfiles

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_gll_impose_val(iregion_code,r,theta,phi,ispec,i,j,k, &
                                  vpv,vph,vsv,vsh,rho,eta_aniso, &
                                  c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                  c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

! sets model parameters for specified GLL point (i,j,k,ispec)

  use constants, only: myrank,IREGION_CRUST_MANTLE,IREGION_OUTER_CORE,IREGION_INNER_CORE

  use shared_parameters, only: ANISOTROPIC_3D_MANTLE

  use model_gll_par

  implicit none

  integer,intent(in) :: iregion_code
  double precision,intent(in) :: r,theta,phi
  integer,intent(in) :: ispec,i,j,k

  double precision,intent(inout) :: vpv,vph,vsv,vsh,rho,eta_aniso
  double precision,intent(inout) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                    c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  ! local parameters
  double precision :: vp,vs
  double precision :: Gc_prime,Gs_prime,mu0

  ! only valid for crust/mantle region at the moment...
  select case(iregion_code)
  case (IREGION_CRUST_MANTLE)
    ! crust/mantle
    ! safety check
    if (ispec > MGLL_V%nspec) then
      call exit_MPI(myrank,'model GLL: ispec too big')
    endif

    ! initializes
    Gc_prime = 0.d0
    Gs_prime = 0.d0

    ! gets model values
    select case(MGLL_TYPE)
    case (1)
      ! isotropic model
      ! takes stored GLL values from file
      ! ( note that these values are non-dimensionalized)
      vp = dble( MGLL_V%vp_new(i,j,k,ispec) )
      vs = dble( MGLL_V%vs_new(i,j,k,ispec) )
      rho = dble( MGLL_V%rho_new(i,j,k,ispec) )

      ! isotropic model
      vpv = vp
      vph = vp
      vsv = vs
      vsh = vs
      rho = rho
      eta_aniso = 1.0d0

    case (2)
      ! transverse model parameters
      ! takes stored GLL values from file
      vph = dble( MGLL_V%vph_new(i,j,k,ispec) )
      vpv = dble( MGLL_V%vpv_new(i,j,k,ispec) )
      vsh = dble( MGLL_V%vsh_new(i,j,k,ispec) )
      vsv = dble( MGLL_V%vsv_new(i,j,k,ispec) )
      rho = dble( MGLL_V%rho_new(i,j,k,ispec) )
      eta_aniso = dble( MGLL_V%eta_new(i,j,k,ispec) )

    case (3)
      ! azimuthal model
      vph = dble( MGLL_V%vph_new(i,j,k,ispec) )
      vpv = dble( MGLL_V%vpv_new(i,j,k,ispec) )
      vsh = dble( MGLL_V%vsh_new(i,j,k,ispec) )
      vsv = dble( MGLL_V%vsv_new(i,j,k,ispec) )
      rho = dble( MGLL_V%rho_new(i,j,k,ispec) )
      eta_aniso = dble( MGLL_V%eta_new(i,j,k,ispec) )
      Gc_prime = dble( MGLL_V%Gc_prime_new(i,j,k,ispec) )
      Gs_prime = dble( MGLL_V%Gs_prime_new(i,j,k,ispec) )
      mu0 = dble( MGLL_V%mu0_new(i,j,k,ispec) )
    case default
      stop 'Invalid MGLL_TYPE, imposing val not implemented yet'
    end select

    ! converts to cij parameters
    if (ANISOTROPIC_3D_MANTLE) then
        ! parameters need to be converted to cijkl
        call model_gll_build_cij(r,theta,phi, &
                                 vph,vpv,vsh,vsv,rho,eta_aniso,Gc_prime,Gs_prime,mu0, &
                                 c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                 c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)
    endif ! ANISOTROPIC_3D_MANTLE

  case (IREGION_OUTER_CORE)
    ! outer core, no GLL model implemented yet
    continue

  case (IREGION_INNER_CORE)
    ! inner core
    ! checks if anything to do
    if (MGLL_V_IC%nspec == 0) return

    ! isotropic model
    ! takes stored GLL values from file
    ! ( note that these values are non-dimensionalized)
    vp = dble( MGLL_V_IC%vp_new(i,j,k,ispec) )
    vs = dble( MGLL_V_IC%vs_new(i,j,k,ispec) )
    rho = dble( MGLL_V_IC%rho_new(i,j,k,ispec) )

    ! isotropic model
    vpv = vp
    vph = vp
    vsv = vs
    vsh = vs
    rho = rho
    eta_aniso = 1.0d0

  end select

  end subroutine model_gll_impose_val

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_gll_build_cij(r,theta,phi, &
                                 vph,vpv,vsh,vsv,rho,eta_aniso,Gc_prime,Gs_prime,mu0, &
                                 c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                 c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  use constants, only: ZERO
  use model_gll_par

  implicit none
  double precision,intent(in) :: r,theta,phi
  double precision,intent(in) :: vpv,vph,vsv,vsh,rho,eta_aniso,Gc_prime,Gs_prime,mu0

  double precision,intent(inout) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                    c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  ! local parameters
  double precision :: A,C,L,N,F,Gc,Gs
  double precision :: Jc,Js,Kc,Ks,Mc,Ms,Bc,Bs,Hc,Hs,Dc,Ds,Ec,Es
  double precision :: r_dummy ! to avoid compilation warning

  ! initializes
  A = ZERO
  C = ZERO
  L = ZERO
  N = ZERO
  F = ZERO

  Gc = ZERO
  Gs = ZERO
  Jc = ZERO
  Js = ZERO
  Kc = ZERO
  Ks = ZERO
  Mc = ZERO
  Ms = ZERO
  Bc = ZERO
  Bs = ZERO
  Hc = ZERO
  Hs = ZERO
  Dc = ZERO
  Ds = ZERO
  Ec = ZERO
  Es = ZERO

  ! to avoid compilation warning
  r_dummy = r

  ! local position (d_ij given in radial direction)

! Ebru: No need to scale elastic tensor as wavespeeds are already scaled before
!       the construction of the tensor.
!  d11 = d11!/scale_GPa
! ..

  select case (MGLL_TYPE)
  case (1)
    ! isotropic model
    !
    ! vpv == vp and vsv == vs and eta == 1
    ! no symmetry axis, no need for rotation (following cij should be rotation invariant)
    !
    ! if equivalent with an isotropic tensor, cij can be set with
    ! Lame parameters: mu = rho * vs**2
    !                  lambda = rho * (vp**2 - 2 vs**2) = rho * vp**2 - 2 mu
    !
    !            then: C11 = C22 = C33 = lambda + 2mu
    !                  C12 = C13 = C23 = lambda
    !                  C44 = C55 = C66 = mu
    c11 = rho * vpv*vpv                  ! lambda + 2 mu
    c12 = rho * (vpv*vpv - 2.d0*vsv*vsv) ! lambda
    c13 = c12                            ! lambda
    c14 = 0.d0
    c15 = 0.d0
    c16 = 0.d0
    c22 = c11                            ! lambda + 2 mu
    c23 = c12                            ! lambda
    c24 = 0.d0
    c25 = 0.d0
    c26 = 0.d0
    c33 = c11                            ! lambda + 2 mu
    c34 = 0.d0
    c35 = 0.d0
    c36 = 0.d0
    c44 = rho * vsv*vsv                  ! mu
    c45 = 0.d0
    c46 = 0.d0
    c55 = c44                            ! mu
    c56 = 0.d0
    c66 = c44                            ! mu

  case (2)
    ! anisotropy based on tiso parameters
    !
    !! cij equivalent with an transversely isotropic elastic tensor
    !
    ! Anderson & Dziewonski, 1982: "Upper mantle anisotropy: evidence from free oscillations", GJR
    ! A = rho * vph**2
    ! C = rho * vpv**2
    ! N = rho * vsh**2
    ! L = rho * vsv**2
    ! F = eta * (A - 2*L)
    !
    ! and therefore (assuming radial axis symmetry, see e.g. Stein & Wysession, chapter 3.6.2)
    ! C11 = A = rho * vph**2
    ! C33 = C = rho * vpv**2
    ! C44 = L = rho * vsv**2
    ! C13 = F = eta * (A - 2*L)
    ! C12 = C11 - 2 C66 = A - 2*N = rho * (vph**2 - 2 * vsh**2)
    ! C22 = C11 = A
    ! C23 = C13 = F
    ! C55 = C44 = L
    ! C66 = N = rho * vsh**2 = (C11-C12)/2
    !
    ! Love parameterization
    A = rho * vph**2
    C = rho * vpv**2
    L = rho * vsv**2
    N = rho * vsh**2
    F = eta_aniso * (A - 2.d0 * L)

    ! local (radial) coordinate system to global SPECFEM reference
    call rotate_tensor_Love_to_global(theta,phi, &
                                      A,C,N,L,F, &
                                      c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                      c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  case (3)
    ! azimuthal anisotropy

    !for normalization of Gc & Gs similar to the kernels:
    !
    !rho_p = ZERO
    !vpv_p = ZERO
    !vph_p = ZERO
    !vsv_p = ZERO
    !vsh_p = ZERO
    !eta_aniso_p = ZERO
    !vp_p = ZERO
    !vs_p = ZERO
    !
    !call model_prem_iso(myrank,r_prem,rho,drhodr,vp,vs,Qkappa,Qmu,idoubling,CRUSTAL,.false.)
    !
    !if (REFERENCE_1D_MODEL == REFERENCE_MODEL_1DREF) then
    !   call model_1dref_broadcast(CRUSTAL)
    !   call model_1dref(r_prem,rho_p,vpv_p,vph_p,vsv_p,vsh_p,eta_aniso_p,Qkappa_p,Qmu_p,iregion_code,CRUSTAL)
    !   ! this case here is only executed for 1D_ref_iso
    !   ! calculates isotropic values
    !   vp_p = sqrt(((8.d0+4.d0*eta_aniso_p)*vph_p*vph_p + 3.d0*vpv_p*vpv_p &
    !                + (8.d0 - 8.d0*eta_aniso_p)*vsv_p*vsv_p)/15.d0)
    !   vs_p = sqrt(((1.d0-2.d0*eta_aniso_p)*vph_p*vph_p + vpv_p*vpv_p &
    !                + 5.d0*vsh_p*vsh_p + (6.d0+4.d0*eta_aniso_p)*vsv_p*vsv_p)/15.d0)
    !endif
    !
    !Gc = Gc_prime * (rho_p*vs_p**2)
    !Gs = Gs_prime * (rho_p*vs_p**2)

    ! note: for Gc_prime and Gs_prime, the scaling rho * vs0**2 is equal to the isotropic shear moduli mu0 = rho * vs0**2
    !       The mu0 value could in principle be choosen arbitrarily, as long as it is coherent with the scaling
    !       for the kernel values. Its purpose is to non-dimensionalize Gc and scale the kernel contributions
    !       with respect to the other ones and balance the contributions.
    !
    ! test: choosing PREM crustal values: rho=2.6 g/cm3, vp=5.8 km/s, vs=3.2 km/s -> mu0 = 26.624 GPa
    !Gc = Gc_prime * 26.6/scale_GPa
    !Gs = Gs_prime * 26.6/scale_GPa

    Gc = Gc_prime * mu0
    Gs = Gs_prime * mu0

    !daniel debug:
    ! test: ignore Gc,Gs contributions
    !Gc = 0.d0
    !Gs = 0.d0

    ! Love parameterization
    A = rho * vph**2    !rhovphsq = A  !!! that is A
    C = rho * vpv**2    !rhovpvsq = C  !!! that is C
    L = rho * vsv**2    !rhovsvsq = L  !!! that is L
    N = rho * vsh**2    !rhovshsq = N  !!! that is N
    F = eta_aniso * (A - 2.d0 * L)

    ! local (azimuthal) coordinate system to global SPECFEM reference
    call rotate_tensor_azimuthal_to_global(theta,phi, &
                                           A,C,N,L,F, &
                                           Gc,Gs, &
                                           c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                           c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)


  case default
    ! saftey stop
    stop 'Invalid MGLL_TYPE, conversion to elastic tensor not implemented yet'

    !! for full parameterization, could call something like:
    !   call rotate_tensor_aniso_to_global(theta,phi &
    !                                      A,C,N,L,F, &
    !                                      Gc,Gs, &
    !                                      Jc,Js,Kc,Ks,Mc,Ms,Bc,Bs,Hc,Hs,Dc,Ds,Ec,Es, &
    !                                      c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
    !                                      c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)
  end select

  ! debug
  !if (myrank==0) print *,"c11",c11

  end subroutine model_gll_build_cij
