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
    double precision :: scale_velocity,scale_density
    ! isotropic model
    real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: vs_new,vp_new,rho_new
    ! transverse isotropic model
    real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: vsv_new,vpv_new, &
      vsh_new,vph_new,eta_new
    ! number of elements (crust/mantle elements)
    integer :: nspec
    integer :: dummy_pad ! padding 4 bytes to align the structure
  end type model_gll_variables
  type (model_gll_variables) :: MGLL_V

  end module model_gll_par

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_gll_broadcast()

! standard routine to setup model

  use constants

  use shared_parameters, only: TRANSVERSE_ISOTROPY,NSPEC_REGIONS,ADIOS_FOR_MODELS,NPROCTOT,NCHUNKS

  use model_gll_par

  implicit none

  ! local parameters
  double precision :: scaleval
  real(kind=CUSTOM_REAL) :: minvalue,maxvalue,min_all,max_all
  integer :: ier

  ! sets number of elements (crust/mantle region)
  MGLL_V%nspec = NSPEC_REGIONS(IREGION_CRUST_MANTLE)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'broadcast model: GLL'
    write(IMAIN,*) 'setup:'
    write(IMAIN,*) '  NCHUNKS           : ',NCHUNKS
    write(IMAIN,*) '  NPROC total       : ',NPROCTOT
    write(IMAIN,*) '  NSPEC             : ',MGLL_V%nspec
    write(IMAIN,*) '  NGLLX/NGLLY/NGLLZ : ',NGLLX,NGLLY,NGLLZ
    call flush_IMAIN()
  endif

  ! allocates arrays
  ! differs for isotropic model or transverse isotropic models
  if (.not. TRANSVERSE_ISOTROPY) then
    ! isotropic model
    allocate( MGLL_V%vp_new(NGLLX,NGLLY,NGLLZ,MGLL_V%nspec), &
              MGLL_V%vs_new(NGLLX,NGLLY,NGLLZ,MGLL_V%nspec), stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating vp_new,.. arrays')
    MGLL_V%vp_new(:,:,:,:) = 0.0_CUSTOM_REAL
    MGLL_V%vs_new(:,:,:,:) = 0.0_CUSTOM_REAL
  else
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
  endif
  allocate( MGLL_V%rho_new(NGLLX,NGLLY,NGLLZ,MGLL_V%nspec), stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating rho_new,.. arrays')
  MGLL_V%rho_new(:,:,:,:) = 0.0_CUSTOM_REAL

  ! reads in model files for each process
  if (ADIOS_FOR_MODELS) then
    call read_gll_model_adios(myrank)
  else
    call read_gll_model(myrank)
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
      call flush_IMAIN()
    endif
    ! Vp
    maxvalue = maxval( MGLL_V%vp_new )
    minvalue = minval( MGLL_V%vp_new )
    call max_all_cr(maxvalue, max_all)
    call min_all_cr(minvalue, min_all)
    if (myrank == 0) then
      write(IMAIN,*) '  vp new min/max: ',min_all,max_all
      call flush_IMAIN()
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
      call flush_IMAIN()
    endif

    ! Vsv
    maxvalue = maxval( MGLL_V%vsv_new )
    minvalue = minval( MGLL_V%vsv_new )
    call max_all_cr(maxvalue, max_all)
    call min_all_cr(minvalue, min_all)
    if (myrank == 0) then
      write(IMAIN,*) '  vsv new min/max: ',min_all,max_all
      call flush_IMAIN()
    endif
    ! Vsh
    maxvalue = maxval( MGLL_V%vsh_new )
    minvalue = minval( MGLL_V%vsh_new )
    call max_all_cr(maxvalue, max_all)
    call min_all_cr(minvalue, min_all)
    if (myrank == 0) then
      write(IMAIN,*) '  vsh new min/max: ',min_all,max_all
      call flush_IMAIN()
    endif
    ! Vpv
    maxvalue = maxval( MGLL_V%vpv_new )
    minvalue = minval( MGLL_V%vpv_new )
    call max_all_cr(maxvalue, max_all)
    call min_all_cr(minvalue, min_all)
    if (myrank == 0) then
      write(IMAIN,*) '  vpv new min/max: ',min_all,max_all
      call flush_IMAIN()
    endif
    ! Vph
    maxvalue = maxval( MGLL_V%vph_new )
    minvalue = minval( MGLL_V%vph_new )
    call max_all_cr(maxvalue, max_all)
    call min_all_cr(minvalue, min_all)
    if (myrank == 0) then
      write(IMAIN,*) '  vph new min/max: ',min_all,max_all
      call flush_IMAIN()
    endif
    ! density
    maxvalue = maxval( MGLL_V%rho_new )
    minvalue = minval( MGLL_V%rho_new )
    call max_all_cr(maxvalue, max_all)
    call min_all_cr(minvalue, min_all)
    if (myrank == 0) then
      write(IMAIN,*) '  rho new min/max: ',min_all,max_all
      call flush_IMAIN()
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
  call synchronize_all()

  end subroutine model_gll_broadcast

!
!-------------------------------------------------------------------------------------------------
!


  subroutine read_gll_model(rank)

  use constants
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
    write(IMAIN,*)'  reading done'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine read_gll_model


!
!-------------------------------------------------------------------------------------------------
!


  subroutine read_gll_modelfiles(rank)

  use constants, only: MAX_STRING_LEN,IMAIN,IIN,PATHNAME_GLL_modeldir
  use shared_parameters, only: TRANSVERSE_ISOTROPY

  use model_gll_par

  implicit none

  integer,intent(in) :: rank

  ! local parameters
  integer :: ier
  character(len=MAX_STRING_LEN) :: prname

  ! for example: DATA/GLL/proc000000_reg1_rho.bin
  !
  ! root name
  write(prname,'(a,i6.6,a)') PATHNAME_GLL_modeldir(1:len_trim(PATHNAME_GLL_modeldir))//'proc',rank,'_reg1_'

  ! reads in model for each partition
  if (.not. TRANSVERSE_ISOTROPY) then
    if (rank == 0) then
      write(IMAIN,*)'  reads isotropic model values: vp,vs,rho'
      call flush_IMAIN()
    endif

    ! isotropic model
    ! vp mesh
    open(unit=IIN,file=prname(1:len_trim(prname))//'vp.bin', &
          status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      write(IMAIN,*) 'Error opening: ',prname(1:len_trim(prname))//'vp.bin'
      call exit_MPI(rank,'Error model GLL')
    endif
    read(IIN) MGLL_V%vp_new(:,:,:,1:MGLL_V%nspec)
    close(IIN)

    ! vs mesh
    open(unit=IIN,file=prname(1:len_trim(prname))//'vs.bin', &
         status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',prname(1:len_trim(prname))//'vs.bin'
      call exit_MPI(rank,'Error model GLL')
    endif
    read(IIN) MGLL_V%vs_new(:,:,:,1:MGLL_V%nspec)
    close(IIN)

  else
    if (rank == 0) then
      write(IMAIN,*)'  reads transversely isotropic model values: vpv,vph,vsv,vsh,eta,rho'
      call flush_IMAIN()
    endif

    ! transverse isotropic model
    ! vp mesh
    open(unit=IIN,file=prname(1:len_trim(prname))//'vpv.bin', &
          status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      write(IMAIN,*) 'Error opening: ',prname(1:len_trim(prname))//'vpv.bin'
      call exit_MPI(rank,'Error model GLL')
    endif
    read(IIN) MGLL_V%vpv_new(:,:,:,1:MGLL_V%nspec)
    close(IIN)

    open(unit=IIN,file=prname(1:len_trim(prname))//'vph.bin', &
          status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      write(IMAIN,*) 'Error opening: ',prname(1:len_trim(prname))//'vph.bin'
      call exit_MPI(rank,'Error model GLL')
    endif
    read(IIN) MGLL_V%vph_new(:,:,:,1:MGLL_V%nspec)
    close(IIN)

    ! vs mesh
    open(unit=IIN,file=prname(1:len_trim(prname))//'vsv.bin', &
         status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',prname(1:len_trim(prname))//'vsv.bin'
      call exit_MPI(rank,'Error model GLL')
    endif
    read(IIN) MGLL_V%vsv_new(:,:,:,1:MGLL_V%nspec)
    close(IIN)

    open(unit=IIN,file=prname(1:len_trim(prname))//'vsh.bin', &
         status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',prname(1:len_trim(prname))//'vsh.bin'
      call exit_MPI(rank,'Error model GLL')
    endif
    read(IIN) MGLL_V%vsh_new(:,:,:,1:MGLL_V%nspec)
    close(IIN)

    ! eta mesh
    open(unit=IIN,file=prname(1:len_trim(prname))//'eta.bin', &
         status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening: ',prname(1:len_trim(prname))//'eta.bin'
      call exit_MPI(rank,'Error model GLL')
    endif
    read(IIN) MGLL_V%eta_new(:,:,:,1:MGLL_V%nspec)
    close(IIN)

  endif

  ! rho mesh
  open(unit=IIN,file=prname(1:len_trim(prname))//'rho.bin', &
       status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',prname(1:len_trim(prname))//'rho.bin'
    call exit_MPI(rank,'Error model GLL')
  endif
  read(IIN) MGLL_V%rho_new(:,:,:,1:MGLL_V%nspec)
  close(IIN)

  end subroutine read_gll_modelfiles

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_gll_impose_val(vpv,vph,vsv,vsh,rho,dvp,eta_aniso,ispec,i,j,k)

! sets model parameters for specified GLL point (i,j,k,ispec)

  use constants, only: myrank

  use shared_parameters, only: TRANSVERSE_ISOTROPY

  use model_gll_par

  implicit none

  double precision,intent(inout) :: vpv,vph,vsv,vsh,rho,dvp,eta_aniso
  integer,intent(in) :: ispec,i,j,k

  ! local parameters
  double precision :: vp,vs

  ! isotropic model
  if (.not. TRANSVERSE_ISOTROPY) then

    !check
    if (ispec > MGLL_V%nspec) then
      call exit_MPI(myrank,'model GLL: ispec too big')
    endif

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

  ! transverse isotropic model
  else

    !check
    if (ispec > MGLL_V%nspec) then
      call exit_MPI(myrank,'model GLL: ispec too big')
    endif

    ! takes stored GLL values from file
    vph = dble( MGLL_V%vph_new(i,j,k,ispec) )
    vpv = dble( MGLL_V%vpv_new(i,j,k,ispec) )
    vsh = dble( MGLL_V%vsh_new(i,j,k,ispec) )
    vsv = dble( MGLL_V%vsv_new(i,j,k,ispec) )
    rho = dble( MGLL_V%rho_new(i,j,k,ispec) )
    eta_aniso = dble( MGLL_V%eta_new(i,j,k,ispec) )

  endif

  ! no mantle vp perturbation
  dvp = 0.0d0

  end subroutine model_gll_impose_val

