program test_models

  use constants, only: DEGREES_TO_RADIANS,MAX_STRING_LEN,myrank
  use meshfem_par
  use manager_adios

  implicit none

  ! local parameters
  integer :: imodel

  integer, parameter :: num_model_names = 56
  character(len=MAX_STRING_LEN),parameter :: model_names(num_model_names) = (/character(len=MAX_STRING_LEN) :: &
    ! 1D models with real structure
    "1D_isotropic_prem", &
    "1D_transversely_isotropic_prem", &
    "1D_isotropic_prem2", &
    "1D_transversely_isotropic_prem2", &
    "1D_iasp91", &
    "1D_1066a", &
    "1D_ak135f_no_mud", &
    "1D_ref", &
    "1D_ref_iso", &
    "1D_jp3d", &
    "1D_sea99", &
    "1D_CCREM", &
    "Ishii", &
    ! 1D models with only one fictitious crustal layer
    "1D_isotropic_prem_onecrust", &
    "1D_transversely_isotropic_prem_onecrust", &
    "1D_isotropic_prem2_onecrust", &
    "1D_transversely_isotropic_prem2_onecrust", &
    "1D_iasp91_onecrust", &
    "1D_1066a_onecrust", &
    "1D_ak135f_no_mud_onecrust", &
    ! fully 3D models
    "transversely_isotropic_prem_plus_3D_crust_2.0", &
    "transversely_isotropic_prem2_plus_3D_crust_2.0", &
    "3D_anisotropic", &
    "3D_attenuation", &
    "s20rts", &
    "s40rts", &
    "s362ani", &
    "s362iso", &
    "s362wmani", &
    "s362ani_prem", &
    "s362ani_3DQ", &
    "s362iso_3DQ", &
    "s29ea", &
    "sea99_jp3d1994", &
    "sea99", &
    "jp3d1994", &
    "sgloberani_aniso", &
    "sgloberani_iso", &
    "SPiRaL", &
    "GLAD_bkmns", &
    "gapp2", &
    ! 3D crustal models append: **mantle_model_name** + _ + crust1.0, crust2.0, EPcrust, EuCRUST, crustmaps, crustSH
    "s20rts_crust1.0", &
    "s20rts_crust2.0", &
    "s20rts_EPcrust", &
    "s20rts_EuCRUST", &
    "s20rts_crustmaps", &
    "s20rts_crustSH", &
    "s20rts_crustSPiRaL", &
    ! 3D models with 1D crust: append "_1Dcrust" to the 3D model name
    "s20rts_1Dcrust", &
    "s362ani_1Dcrust", &
    ! generic models
    "PPM", &
    "full_sh", &
    "heterogen", &
    ! Mars
    "1D_Sohl", &
    "1D_case65TAY", &
    ! Moon
    "VPREMOON" &
    /)

  ! not tested: cem_request, cem_accept, cem_gll
  !             would need compilation --with-cem

  ! initialize the MPI communicator and start the NPROCTOT MPI processes.
  call init_mpi()
  call world_rank(myrank)

  if (myrank == 0) print *,'program: test_models'

  ! first try read the Par_file
  if (myrank == 0) then
    ! reads the parameter file
    call read_parameter_file()
  endif

  ! check different models
  if (myrank == 0) print *,'testing model setup...'
  do imodel = 1,num_model_names
    ! reads the parameter file and computes additional parameters
    call read_parameter_file()

    ! re-set model name
    MODEL = model_names(imodel)

    ! specific setting for models: Mars & Moon models cannot have OCEANS
    if (trim(MODEL) == "1D_Sohl" &
        .or. trim(MODEL) == "1D_case65TAY" &
        .or. trim(MODEL) == "VPREMOON") then
      OCEANS = .false.
    endif

    ! tests model setups
    call test_initialize_models()
  enddo

  ! tests further calls from initialize_mesher.f90

  ! compute rotation matrix from Euler angles
  ANGULAR_WIDTH_XI_RAD = ANGULAR_WIDTH_XI_IN_DEGREES * DEGREES_TO_RADIANS
  ANGULAR_WIDTH_ETA_RAD = ANGULAR_WIDTH_ETA_IN_DEGREES * DEGREES_TO_RADIANS

  if (NCHUNKS /= 6) call euler_angles(rotation_matrix,CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH)

  ! ADIOS
  if (ADIOS_ENABLED) then
    call initialize_adios()
  endif

  ! gravity integrals
  if (GRAVITY_INTEGRALS) then
    call gravity_initialize_integrals()
  endif

  ! OpenMP
  if (myrank == 0) print *,'testing OpenMP setup...'
  call init_openmp()

  ! done
  if (myrank == 0) print *,'test_models done successfully'

  ! stop all the MPI processes, and exit
  call finalize_mpi()

end program test_models

!
!-------------------------------------------------------------------------------
!

  subroutine test_initialize_models()

  use meshfem_par
  use meshfem_models_par

  implicit none

  ! user output
  if (myrank == 0) then
    print *,'MODEL: ',trim(MODEL)
  endif

  if (myrank == 0) then
    ! sets compute parameters accordingly
    call rcp_set_compute_parameters()

    ! sets mesh parameters
    call rcp_set_mesh_parameters()
  endif

  ! broadcast parameters read from main process to all processes
  call broadcast_computed_parameters()

  ! synchronizes processes
  call synchronize_all()

  end subroutine test_initialize_models

