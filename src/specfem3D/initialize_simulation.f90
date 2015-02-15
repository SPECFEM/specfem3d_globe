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

  subroutine initialize_simulation()

  use specfem_par
  use specfem_par_movie

  implicit none

  ! local parameters
  integer :: sizeprocs
  integer :: ier
  character(len=MAX_STRING_LEN) :: dummystring
  character(len=MAX_STRING_LEN) :: path_to_add

  ! sizeprocs returns number of processes started (should be equal to NPROCTOT).
  ! myrank is the rank of each process, between 0 and sizeprocs-1.
  ! as usual in MPI, process 0 is in charge of coordinating everything
  ! and also takes care of the main output
  call world_size(sizeprocs)
  call world_rank(myrank)

  ! open main output file, only written to by process 0
  if (myrank == 0) then
    if (IMAIN /= ISTANDARD_OUTPUT) then
      open(unit=IMAIN,file=trim(OUTPUT_FILES)//'/output_solver.txt',status='unknown',action='write',iostat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error opening file output_solver.txt for writing output info')
    endif

    write(IMAIN,*)
    write(IMAIN,*) '******************************'
    write(IMAIN,*) '**** Specfem3D MPI Solver ****'
    write(IMAIN,*) '******************************'
    write(IMAIN,*)
    write(IMAIN,*)
    call flush_IMAIN()
  endif

!! DK DK for Roland_Sylvain
  if (ROLAND_SYLVAIN) call exit_MPI(myrank,'no need to run the solver to compute Roland_Sylvain integrals, only the mesher')

  if (myrank == 0) then
    ! read the parameter file and compute additional parameters
    call read_compute_parameters()

!! DK DK make sure NSTEP is a multiple of NT_DUMP_ATTENUATION
!! DK DK we cannot move this to inside read_compute_parameters because when read_compute_parameters
!! DK DK is called from the beginning of create_header_file then the value of NT_DUMP_ATTENUATION is unknown
    if (UNDO_ATTENUATION .and. mod(NSTEP,NT_DUMP_ATTENUATION) /= 0) then
      NSTEP = (NSTEP/NT_DUMP_ATTENUATION + 1)*NT_DUMP_ATTENUATION
      ! subsets used to save seismograms must not be larger than the whole time series, otherwise we waste memory
      if (NTSTEP_BETWEEN_OUTPUT_SEISMOS > NSTEP) NTSTEP_BETWEEN_OUTPUT_SEISMOS = NSTEP
    endif
  endif

  ! broadcast parameters read from master to all processes
  call broadcast_computed_parameters(myrank)

  ! check that the code is running with the requested nb of processes
  if (sizeprocs /= NPROCTOT) then
    print*,'Error: rank ',myrank,' - wrong number of MPI processes',sizeprocs,NPROCTOT
    call exit_MPI(myrank,'wrong number of MPI processes in the initialization of SPECFEM')
  endif

  ! synchronizes processes
  call synchronize_all()

  if (myrank == 0) then

    if (FIX_UNDERFLOW_PROBLEM) write(IMAIN,*) 'Fixing slow underflow trapping problem using small initial field'

    write(IMAIN,*)
    write(IMAIN,*) 'There are ',sizeprocs,' MPI processes'
    write(IMAIN,*) 'Processes are numbered from 0 to ',sizeprocs-1
    write(IMAIN,*)

    write(IMAIN,*) 'There are ',NEX_XI,' elements along xi in each chunk'
    write(IMAIN,*) 'There are ',NEX_ETA,' elements along eta in each chunk'
    write(IMAIN,*)
    write(IMAIN,*) 'There are ',NPROC_XI,' slices along xi in each chunk'
    write(IMAIN,*) 'There are ',NPROC_ETA,' slices along eta in each chunk'
    write(IMAIN,*) 'There is a total of ',NPROC,' slices in each chunk'
    write(IMAIN,*) 'There are ',NCHUNKS,' chunks'
    write(IMAIN,*) 'There is a total of ',NPROCTOT,' slices in all the chunks'

    write(IMAIN,*)
    write(IMAIN,*) 'NDIM = ',NDIM
    write(IMAIN,*)
    write(IMAIN,*) 'NGLLX = ',NGLLX
    write(IMAIN,*) 'NGLLY = ',NGLLY
    write(IMAIN,*) 'NGLLZ = ',NGLLZ
    write(IMAIN,*)

    ! write information about precision used for floating-point operations
    if (CUSTOM_REAL == SIZE_REAL) then
      write(IMAIN,*) 'using single precision for the calculations'
    else
      write(IMAIN,*) 'using double precision for the calculations'
    endif
    write(IMAIN,*)
    write(IMAIN,*) 'smallest and largest possible floating-point numbers are: ', &
      tiny(1._CUSTOM_REAL),huge(1._CUSTOM_REAL)
    write(IMAIN,*)

    ! model user parameters
    write(IMAIN,*) 'model: ',trim(MODEL)
    if (OCEANS) then
      write(IMAIN,*) '  incorporating the oceans using equivalent load'
    else
      write(IMAIN,*) '  no oceans'
    endif
    if (ELLIPTICITY) then
      write(IMAIN,*) '  incorporating ellipticity'
    else
      write(IMAIN,*) '  no ellipticity'
    endif
    if (TOPOGRAPHY) then
      write(IMAIN,*) '  incorporating surface topography'
    else
      write(IMAIN,*) '  no surface topography'
    endif
    if (GRAVITY) then
      write(IMAIN,*) '  incorporating self-gravitation (Cowling approximation)'
    else
      write(IMAIN,*) '  no self-gravitation'
    endif
    if (ROTATION) then
      write(IMAIN,*) '  incorporating rotation'
    else
      write(IMAIN,*) '  no rotation'
    endif
    if (ATTENUATION) then
      write(IMAIN,*) '  incorporating attenuation using ',N_SLS,' standard linear solids'
      if (ATTENUATION_3D) write(IMAIN,*)'  using 3D attenuation model'
    else
      write(IMAIN,*) '  no attenuation'
    endif
    write(IMAIN,*)

    ! model mesh parameters
    if (ISOTROPIC_3D_MANTLE) then
      write(IMAIN,*) '  incorporating 3-D lateral variations'
    else
      write(IMAIN,*) '  no 3-D lateral variations'
    endif
    if (HETEROGEN_3D_MANTLE) then
      write(IMAIN,*) '  incorporating heterogeneities in the mantle'
    else
      write(IMAIN,*) '  no heterogeneities in the mantle'
    endif
    if (CRUSTAL) then
      write(IMAIN,*) '  incorporating crustal variations'
    else
      write(IMAIN,*) '  no crustal variations'
    endif
    if (ONE_CRUST) then
      write(IMAIN,*) '  using one layer only in PREM crust'
    else
      write(IMAIN,*) '  using unmodified 1D crustal model with two layers'
    endif
    if (TRANSVERSE_ISOTROPY) then
      write(IMAIN,*) '  incorporating transverse isotropy'
    else
      write(IMAIN,*) '  no transverse isotropy'
    endif
    if (ANISOTROPIC_INNER_CORE_VAL) then
      write(IMAIN,*) '  incorporating anisotropic inner core'
    else
      write(IMAIN,*) '  no inner-core anisotropy'
    endif
    if (ANISOTROPIC_3D_MANTLE_VAL) then
      write(IMAIN,*) '  incorporating anisotropic mantle'
    else
      write(IMAIN,*) '  no general mantle anisotropy'
    endif

    write(IMAIN,*)
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! define strain storage
  ! this cannot be made a constant stored in values_from_mesher.h because it depends on SIMULATION_TYPE
  if (ATTENUATION_VAL .or. SIMULATION_TYPE /= 1 .or. SAVE_FORWARD &
    .or. (MOVIE_VOLUME .and. SIMULATION_TYPE /= 3)) then
    COMPUTE_AND_STORE_STRAIN = .true.
  else
    COMPUTE_AND_STORE_STRAIN = .false.
  endif

  ! checks flags
  call initialize_simulation_check()

  ! synchronizes processes
  call synchronize_all()

  ! counts receiver stations
  if (SIMULATION_TYPE == 1) then
    rec_filename = 'DATA/STATIONS'
  else
    rec_filename = 'DATA/STATIONS_ADJOINT'
  endif
  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
    write(path_to_add,"('run',i4.4,'/')") mygroup + 1
    rec_filename=path_to_add(1:len_trim(path_to_add))//rec_filename(1:len_trim(rec_filename))
  endif
  STATIONS = rec_filename

  ! get total number of receivers
  if (myrank == 0) then
    open(unit=IIN,file=STATIONS,iostat=ier,status='old',action='read')
    nrec = 0
    do while(ier == 0)
      read(IIN,"(a)",iostat=ier) dummystring
      if (ier == 0) then
        ! excludes empty lines
        if (len_trim(dummystring) > 0 ) nrec = nrec + 1
      endif
    enddo
    close(IIN)
  endif

  ! broadcast the information read on the master to the nodes
  call bcast_all_singlei(nrec)

  ! checks number of total receivers
  if (nrec < 1) call exit_MPI(myrank,trim(STATIONS)//': need at least one receiver')

  ! initializes GPU cards
  call initialize_GPU()

  ! initializes VTK window
  if (VTK_MODE) then
    if (myrank == 0 ) call initialize_vtkwindow(GPU_MODE)
  endif

  if (ADIOS_ENABLED .or. OUTPUT_SEISMOS_ASDF) then
    call adios_setup()
  endif
  !if (ADIOS_ENABLED) then
    ! TODO use only one ADIOS group to write simulation parameters
    !      i.e. merge write_solver... write_par_... into
    !      write_specfem3D_globe_adios_header()
    !call write_solver_info_header_ADIOS()
    !call write_specfem_header_adios()
  !endif

  ! synchronizes processes
  call synchronize_all()

  end subroutine initialize_simulation

!
!-------------------------------------------------------------------------------------------------
!

  subroutine initialize_simulation_check()

  use specfem_par
  implicit none

  ! check that the code has been compiled with the right values
  if (NSPEC(IREGION_CRUST_MANTLE) /= NSPEC_CRUST_MANTLE) then
      if (myrank == 0) write(IMAIN,*) 'NSPEC_CRUST_MANTLE:',NSPEC(IREGION_CRUST_MANTLE),NSPEC_CRUST_MANTLE
      write(*,*) 'NSPEC_CRUST_MANTLE:', NSPEC(IREGION_CRUST_MANTLE), NSPEC_CRUST_MANTLE
      call exit_MPI(myrank,'Error in compiled parameters, please recompile solver 1')
  endif
  if (NSPEC(IREGION_OUTER_CORE) /= NSPEC_OUTER_CORE) then
      if (myrank == 0) write(IMAIN,*) 'NSPEC_OUTER_CORE:',NSPEC(IREGION_OUTER_CORE),NSPEC_OUTER_CORE
      write(*,*) 'NSPEC_OUTER_CORE:', NSPEC(IREGION_OUTER_CORE), NSPEC_OUTER_CORE
      call exit_MPI(myrank,'Error in compiled parameters, please recompile solver 2')
  endif
  if (NSPEC(IREGION_INNER_CORE) /= NSPEC_INNER_CORE) then
      if (myrank == 0) write(IMAIN,*) 'NSPEC_INNER_CORE:',NSPEC(IREGION_INNER_CORE),NSPEC_INNER_CORE
      write(*,*) 'NSPEC_INNER_CORE:', NSPEC(IREGION_INNER_CORE), NSPEC_INNER_CORE
      call exit_MPI(myrank,'Error in compiled parameters, please recompile solver 3')
  endif
  if (ATTENUATION_3D .NEQV. ATTENUATION_3D_VAL) then
      if (myrank == 0) write(IMAIN,*) 'ATTENUATION_3D:',ATTENUATION_3D,ATTENUATION_3D_VAL
      write(*,*) 'ATTENUATION_3D:', ATTENUATION_3D, ATTENUATION_3D_VAL
      call exit_MPI(myrank,'Error in compiled parameters ATTENUATION_3D, please recompile solver')
  endif
  if (NCHUNKS /= NCHUNKS_VAL) then
      if (myrank == 0) write(IMAIN,*) 'NCHUNKS:',NCHUNKS,NCHUNKS_VAL
      write(*,*) 'NCHUNKS:', NCHUNKS, NCHUNKS_VAL
      call exit_MPI(myrank,'Error in compiled parameters NCHUNKS, please recompile solver')
  endif
  if (GRAVITY .NEQV. GRAVITY_VAL) then
      if (myrank == 0) write(IMAIN,*) 'GRAVITY:',GRAVITY,GRAVITY_VAL
      write(*,*) 'GRAVITY:', GRAVITY, GRAVITY_VAL
      call exit_MPI(myrank,'Error in compiled parameters GRAVITY, please recompile solver')
  endif
  if (ROTATION .NEQV. ROTATION_VAL) then
      if (myrank == 0) write(IMAIN,*) 'ROTATION:',ROTATION,ROTATION_VAL
      write(*,*) 'ROTATION:', ROTATION, ROTATION_VAL
      call exit_MPI(myrank,'Error in compiled parameters ROTATION, please recompile solver')
  endif
  if (ATTENUATION .NEQV. ATTENUATION_VAL) then
      if (myrank == 0) write(IMAIN,*) 'ATTENUATION:',ATTENUATION,ATTENUATION_VAL
      write(*,*) 'ATTENUATION:', ATTENUATION, ATTENUATION_VAL
      call exit_MPI(myrank,'Error in compiled parameters ATTENUATION, please recompile solver')
  endif
  if (ELLIPTICITY .NEQV. ELLIPTICITY_VAL) then
      if (myrank == 0) write(IMAIN,*) 'ELLIPTICITY:',ELLIPTICITY,ELLIPTICITY_VAL
      write(*,*) 'ELLIPTICITY:', ELLIPTICITY, ELLIPTICITY_VAL
      call exit_MPI(myrank,'Error in compiled parameters ELLIPTICITY, please recompile solver')
  endif
  if (OCEANS .NEQV. OCEANS_VAL) then
      if (myrank == 0) write(IMAIN,*) 'OCEANS:',OCEANS,OCEANS_VAL
      write(*,*) 'OCEANS:', OCEANS, OCEANS_VAL
      call exit_MPI(myrank,'Error in compiled parameters OCEANS, please recompile solver')
  endif
  if (NPROC_XI /= NPROC_XI_VAL) then
      if (myrank == 0) write(IMAIN,*) 'NPROC_XI:',NPROC_XI,NPROC_XI_VAL
      write(*,*) 'NPROC_XI:', NPROC_XI, NPROC_XI_VAL
      call exit_MPI(myrank,'Error in compiled parameters NPROC_XI, please recompile solver')
  endif
  if (NPROC_ETA /= NPROC_ETA_VAL) then
      if (myrank == 0) write(IMAIN,*) 'NPROC_ETA:',NPROC_ETA,NPROC_ETA_VAL
      write(*,*) 'NPROC_ETA:', NPROC_ETA, NPROC_ETA_VAL
      call exit_MPI(myrank,'Error in compiled parameters NPROC_ETA, please recompile solver')
  endif
  if (NPROCTOT /= NPROCTOT_VAL) then
      if (myrank == 0) write(IMAIN,*) 'NPROCTOT:',NPROCTOT,NPROCTOT_VAL
      write(*,*) 'NPROCTOT:', NPROCTOT, NPROCTOT_VAL
      call exit_MPI(myrank,'Error in compiled parameters NPROCTOT, please recompile solver')
  endif
  if (NEX_XI /= NEX_XI_VAL) then
      if (myrank == 0) write(IMAIN,*) 'NEX_XI:',NEX_XI,NEX_XI_VAL
      write(*,*) 'NEX_XI:', NEX_XI, NEX_XI_VAL
      call exit_MPI(myrank,'Error in compiled parameters NEX_XI, please recompile solver')
  endif
  if (NEX_ETA /= NEX_ETA_VAL) then
      if (myrank == 0) write(IMAIN,*) 'NEX_ETA:',NEX_ETA,NEX_ETA_VAL
      write(*,*) 'NEX_ETA:', NEX_ETA, NEX_ETA_VAL
      call exit_MPI(myrank,'Error in compiled parameters NEX_ETA, please recompile solver')
  endif
  if (TRANSVERSE_ISOTROPY .NEQV. TRANSVERSE_ISOTROPY_VAL) then
      if (myrank == 0) write(IMAIN,*) 'TRANSVERSE_ISOTROPY:',TRANSVERSE_ISOTROPY,TRANSVERSE_ISOTROPY_VAL
      write(*,*) 'TRANSVERSE_ISOTROPY:', TRANSVERSE_ISOTROPY, TRANSVERSE_ISOTROPY_VAL
      call exit_MPI(myrank,'Error in compiled parameters, please recompile solver 14')
  endif
  if (ANISOTROPIC_3D_MANTLE .NEQV. ANISOTROPIC_3D_MANTLE_VAL) then
      if (myrank == 0) write(IMAIN,*) 'ANISOTROPIC_3D_MANTLE:',ANISOTROPIC_3D_MANTLE,ANISOTROPIC_3D_MANTLE_VAL
      write(*,*) 'ANISOTROPIC_3D_MANTLE:', ANISOTROPIC_3D_MANTLE, ANISOTROPIC_3D_MANTLE_VAL
      call exit_MPI(myrank,'Error in compiled parameters, please recompile solver 15')
  endif
  if (ANISOTROPIC_INNER_CORE .NEQV. ANISOTROPIC_INNER_CORE_VAL) then
      if (myrank == 0) write(IMAIN,*) 'ANISOTROPIC_INNER_CORE:',ANISOTROPIC_INNER_CORE,ANISOTROPIC_INNER_CORE_VAL
      write(*,*) 'ANISOTROPIC_INNER_CORE:', ANISOTROPIC_INNER_CORE, ANISOTROPIC_INNER_CORE_VAL
      call exit_MPI(myrank,'Error in compiled parameters, please recompile solver 16')
  endif

  ! check simulation parameters
  if (SIMULATION_TYPE /= 1 .and.  SIMULATION_TYPE /= 2 .and. SIMULATION_TYPE /= 3) &
    call exit_MPI(myrank, 'SIMULATION_TYPE can only be 1, 2, or 3')

  ! checks number of sources for adjoint simulations
  ! The limit below is somewhat arbitrary. For pure adjoint simulations (SIMULATION_TYPE == 2),
  ! the code outputs displacement (NT.S00001.BXX.semd,..) and strains (NT.S00001.SEE.semd,..)
  ! as well as source derivative kernels (src_frechet.00001,..) all for each point source.
  ! The naming convention for these files uses (.., i6.6,..), which limits the number of sources to 999999.
  ! If that is still too low, you can increase it further (if so, change all the occurrences of (.., i6.6,..) in the code).
  if (SIMULATION_TYPE /= 1 .and. NSOURCES > 999999) &
    call exit_MPI(myrank,'for adjoint simulations, NSOURCES <= 999999, if you need more change i6.6 in write_seismograms.f90')

  ! attenuation
  if (UNDO_ATTENUATION .and. PARTIAL_PHYS_DISPERSION_ONLY) &
    call exit_MPI(myrank,'cannot have both UNDO_ATTENUATION and PARTIAL_PHYS_DISPERSION_ONLY, please check Par_file...')

  if ((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then

    if (ATTENUATION_VAL) then
      ! checks mimic flag:
      ! attenuation for adjoint simulations must have PARTIAL_PHYS_DISPERSION_ONLY set by xcreate_header_file
      if (.not. EXACT_UNDOING_TO_DISK) then
        if (.not. UNDO_ATTENUATION) then
          if (.not. PARTIAL_PHYS_DISPERSION_ONLY) then
            call exit_MPI(myrank, &
                    'ATTENUATION for adjoint runs or SAVE_FORWARD requires UNDO_ATTENUATION or PARTIAL_PHYS_DISPERSION_ONLY')
          endif
        endif
      endif

      if (PARTIAL_PHYS_DISPERSION_ONLY .NEQV. PARTIAL_PHYS_DISPERSION_ONLY_VAL) then
        if (myrank == 0) write(IMAIN,*) 'PARTIAL_PHYS_DISPERSION_ONLY:',PARTIAL_PHYS_DISPERSION_ONLY, &
                                                                       PARTIAL_PHYS_DISPERSION_ONLY_VAL
        write(*,*) 'PARTIAL_PHYS_DISPERSION_ONLY:', PARTIAL_PHYS_DISPERSION_ONLY, &
                                                    PARTIAL_PHYS_DISPERSION_ONLY_VAL
        call exit_MPI(myrank,'Error in compiled parameters, please recompile solver 17')
      endif

      ! user output
      if (myrank == 0) then
        write(IMAIN,*) 'incorporates ATTENUATION for time-reversed simulation'
        write(IMAIN,*)
        call flush_IMAIN()
      endif
    endif

    ! checks adjoint array dimensions
    if (NSPEC_CRUST_MANTLE_ADJOINT /= NSPEC_CRUST_MANTLE &
      .or. NSPEC_OUTER_CORE_ADJOINT /= NSPEC_OUTER_CORE &
      .or. NSPEC_INNER_CORE_ADJOINT /= NSPEC_INNER_CORE &
      .or. NGLOB_CRUST_MANTLE_ADJOINT /= NGLOB_CRUST_MANTLE &
      .or. NGLOB_OUTER_CORE_ADJOINT /= NGLOB_OUTER_CORE &
      .or. NGLOB_INNER_CORE_ADJOINT /= NGLOB_INNER_CORE) &
      call exit_MPI(myrank, 'improper dimensions of adjoint arrays, please recompile solver 18')
  endif

  ! checks attenuation
  if (ATTENUATION_VAL) then
    if (NSPEC_CRUST_MANTLE_ATTENUATION /= NSPEC_CRUST_MANTLE) &
       call exit_MPI(myrank, 'NSPEC_CRUST_MANTLE_ATTENUATION /= NSPEC_CRUST_MANTLE, exit')
    if (NSPEC_INNER_CORE_ATTENUATION /= NSPEC_INNER_CORE) &
       call exit_MPI(myrank, 'NSPEC_INNER_CORE_ATTENUATION /= NSPEC_INNER_CORE, exit')
  endif

  ! checks strain storage
  if (ATTENUATION_VAL .or. SIMULATION_TYPE /= 1 .or. SAVE_FORWARD &
    .or. (MOVIE_VOLUME .and. SIMULATION_TYPE /= 3)) then
    if (COMPUTE_AND_STORE_STRAIN .neqv. .true. ) &
      call exit_MPI(myrank, 'Error in compiled COMPUTE_AND_STORE_STRAIN parameter, please recompile solver 19')
  else
    if (COMPUTE_AND_STORE_STRAIN .neqv. .false. ) &
      call exit_MPI(myrank, 'Error in compiled COMPUTE_AND_STORE_STRAIN parameter, please recompile solver 20')
  endif

  if (FORCE_VECTORIZATION_VAL) then
    if (ATTENUATION_VAL .and. .not. PARTIAL_PHYS_DISPERSION_ONLY_VAL .and. N_SLS /= 3) &
      call exit_MPI(myrank, &
      'FORCE_VECTORIZATION can only be used with N_SLS == 3 when ATTENUATION .and. .not. PARTIAL_PHYS_DISPERSION_ONLY&
       & because N_SLS is assumed to be equal to 3 when vectorizing compute_element_iso,tiso,aniso')
  endif

  if (SIMULATION_TYPE == 3 .and. (ANISOTROPIC_3D_MANTLE_VAL .or. ANISOTROPIC_INNER_CORE_VAL)) &
     call exit_MPI(myrank, 'anisotropic model is not implemented for kernel simulations yet')

  ! checks model for transverse isotropic kernel computation
  if (SAVE_TRANSVERSE_KL_ONLY) then
    if (ANISOTROPIC_3D_MANTLE_VAL) then
        call exit_mpi(myrank,'Error SAVE_TRANSVERSE_KL_ONLY: Earth model not supported yet')
    endif
    if (SIMULATION_TYPE == 3) then
      if (.not. ANISOTROPIC_KL) then
        call exit_mpi(myrank,'Error SAVE_TRANSVERSE_KL_ONLY: needs anisotropic kernel calculations')
      endif
    endif
  endif

  ! check for GPU runs
  if (GPU_MODE) then
    if (NGLLX /= 5 .or. NGLLY /= 5 .or. NGLLZ /= 5 ) &
      call exit_mpi(myrank,'GPU mode can only be used if NGLLX == NGLLY == NGLLZ == 5')
    if (CUSTOM_REAL /= SIZE_REAL ) &
      call exit_mpi(myrank,'GPU mode runs only with CUSTOM_REAL == SIZE_REAL')
    if (ATTENUATION_VAL) then
      if (N_SLS /= 3 ) &
        call exit_mpi(myrank,'GPU mode does not support N_SLS /= 3 yet')
    endif
  endif

  end subroutine initialize_simulation_check

!
!-------------------------------------------------------------------------------------------------
!

  subroutine initialize_GPU()

! initialization for GPU cards

  use specfem_par
  implicit none
  ! local parameters
  integer :: ngpu_devices,ngpu_devices_min,ngpu_devices_max
  integer :: iproc

  !----------------------------------------------------------------
  ! user test parameters
  !
  ! for hybrid computing: distributes MPI processes to use CPU and GPU
  ! note that a single MPI process on GPU is about >30x faster than on single CPU-core
  !
  ! cray xk7 node: 16-core CPU, 1 K20x GPU card
  !                using 15 processes on single GPU and 1 processes on CPU still slows down the computation
  !                as the GPU slices will be waiting for the CPU one...
  !
  ! turns on/off hybrid CPU-GPU computing
  logical,parameter :: USE_HYBRID_CPU_GPU = .false.
  ! total number of MPI processes run on a single node
  integer, parameter :: TOTAL_PROCESSES_PER_NODE = 16
  ! number of MPI processes run on CPU-cores
  integer, parameter :: PROCESSES_PER_CPU = 1

  !----------------------------------------------------------------

  if (GPU_MODE .and. USE_HYBRID_CPU_GPU) then
    ! distributes processes on GPU and CPU
    if (mod(myrank,TOTAL_PROCESSES_PER_NODE) < PROCESSES_PER_CPU) then
      ! turns of GPU mode for this process
      GPU_MODE = .false.
    else
      !leaves GPU_MODE == .true.
      continue
    endif

    ! user output
    if (myrank == 0 ) print*,'Hybrid CPU-GPU computation:'
    do iproc = 0, NPROCTOT_VAL-1
      if (myrank == iproc) then
        if (myrank < TOTAL_PROCESSES_PER_NODE) then
          print*,'rank ',myrank,' has GPU_MODE set to ',GPU_MODE
        endif
      endif
      call synchronize_all()
    enddo
  endif

  ! initializes number of local gpu devices
  ngpu_devices = 0

  ! GPU_MODE now defined in Par_file
  if (GPU_MODE) then
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) "GPU_MODE Active."
      write(IMAIN,*) "  runtime : ",GPU_RUNTIME
      write(IMAIN,*) "  platform: ",trim(GPU_PLATFORM)
      write(IMAIN,*) "  device  : ",trim(GPU_DEVICE)
      call flush_IMAIN()
    endif

    ! check for GPU runs
    if (NGLLX /= 5 .or. NGLLY /= 5 .or. NGLLZ /= 5 ) &
      stop 'GPU mode can only be used if NGLLX == NGLLY == NGLLZ == 5'
    if (CUSTOM_REAL /= SIZE_REAL) &
      stop 'GPU mode runs only with CUSTOM_REAL == SIZE_REAL'
    if (ATTENUATION_VAL) then
      if (N_SLS /= 3 ) &
        stop 'GPU mode does not support N_SLS /= 3 yet'
    endif

    ! initializes GPU and outputs info to files for all processes
    call initialize_gpu_device(GPU_RUNTIME, GPU_PLATFORM, GPU_DEVICE, myrank,ngpu_devices)
  endif

  ! collects min/max of local devices found for statistics
  call synchronize_all()
  call min_all_i(ngpu_devices,ngpu_devices_min)
  call max_all_i(ngpu_devices,ngpu_devices_max)

  if (GPU_MODE) then
    if (myrank == 0) then
      write(IMAIN,*) "GPU number of devices per node: min =",ngpu_devices_min
      write(IMAIN,*) "                                max =",ngpu_devices_max
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  end subroutine initialize_GPU
