! ADJOINT
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: b_alphaval, b_betaval, b_gammaval

  real(kind=CUSTOM_REAL), dimension(5,N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STR_AND_ATT) :: b_R_memory_crust_mantle
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: b_epsilondev_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: b_eps_trace_over_3_crust_mantle

  real(kind=CUSTOM_REAL), dimension(5,N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_STR_AND_ATT) :: b_R_memory_inner_core
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ADJOINT) :: b_epsilondev_inner_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ADJOINT) :: b_eps_trace_over_3_inner_core

  double precision, dimension(:,:), allocatable :: b_buffer_slices
  double precision, dimension(:,:,:), allocatable :: b_buffer_all_cube_from_slices

! ADJOINT
  real(kind=CUSTOM_REAL) b_deltat,b_deltatover2,b_deltatsqover2
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE_ADJOINT) :: &
    b_displ_crust_mantle,b_veloc_crust_mantle,b_accel_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE_ADJOINT) :: &
    b_displ_outer_core,b_veloc_outer_core,b_accel_outer_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE_ADJOINT) :: &
    b_displ_inner_core,b_veloc_inner_core,b_accel_inner_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ADJOINT) :: b_div_displ_outer_core

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_OUTER_CORE_ADJOINT) :: b_vector_displ_outer_core

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_OUTER_CORE_ADJOINT) :: b_vector_displ_outer_core !ZN

  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_buffer_send_faces,b_buffer_received_faces

  integer :: b_iphase,b_iphase_CC,b_icall

! buffers for send and receive between corners of the chunks
  real(kind=CUSTOM_REAL), dimension(NGLOB1D_RADIAL_CM) :: b_buffer_send_chunkcorn_scalar,b_buffer_recv_chunkcorn_scalar

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB1D_RADIAL_CM + NGLOB1D_RADIAL_IC) :: &
     b_buffer_send_chunkcorn_vector,b_buffer_recv_chunkcorn_vector

!ADJOINT
  real(kind=CUSTOM_REAL) b_two_omega_earth
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ROT_ADJOINT) :: &
    b_A_array_rotation,b_B_array_rotation
