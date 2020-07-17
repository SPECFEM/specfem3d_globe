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
!
!  This portion of the SPECFEM3D Code was written by:
!  Brian Savage while at
!     California Institute of Technology
!     Department of Terrestrial Magnetism / Carnegie Institute of Washington
!     University of Rhode Island
!
!   It is based upon formulation in the following references:
!
!   Dahlen and Tromp, 1998
!      Theoretical Global Seismology
!
!   Liu et al. 1976
!      Velocity dispersion due to anelasticity: implications for seismology and mantle composition
!      Geophys, J. R. asts. Soc, Vol 47, pp. 41-58
!
!   The methodology can be found in Brian Savage, Dimitri Komatitsch and Jeroen Tromp,
!   Effects of 3D attenuation on seismic wave amplitude and phase measurements, Bulletin of the Seismological Society of America,
!   vol. 100(3), p. 1241-1251, doi: 10.1785/0120090263 (2010).
!
! @ARTICLE{SaKoTr10,
!   author = {Brian Savage and Dimitri Komatitsch and Jeroen Tromp},
!   title = {Effects of {3D} attenuation on seismic wave amplitude and phase measurements},
!   journal = {Bull. Seismol. Soc. Am.},
!   year = {2010},
!   volume = {100},
!   pages = {1241-1251},
!   number = {3},
!   doi = {10.1785/0120090263}}
!
!--------------------------------------------------------------------------------------------------

  module model_attenuation_par

  implicit none

  ! model_attenuation_storage_var
  type model_attenuation_storage_var
    sequence
    double precision, dimension(:,:), allocatable :: tau_e_storage
    double precision, dimension(:), allocatable   :: Qmu_storage
    integer :: Q_resolution
    integer :: idummy         ! padding 4 bytes to align the structure
  end type model_attenuation_storage_var
  type (model_attenuation_storage_var) :: AM_S
  ! model_attenuation_storage_var

  ! attenuation_simplex_variables
  type attenuation_simplex_variables
    sequence
    double precision :: Q  ! Q     = Desired Value of Attenuation or Q
    double precision :: iQ ! iQ    = 1/Q
    double precision, dimension(:), allocatable ::  f
    ! f = Frequencies at which to evaluate the solution
    double precision, dimension(:), allocatable :: tau_s
    ! tau_s = Tau_sigma defined by the frequency range and
    !             number of standard linear solids
    integer :: nf          ! nf    = Number of Frequencies
    integer :: nsls        ! nsls  = Number of Standard Linear Solids
  end type attenuation_simplex_variables
  type(attenuation_simplex_variables) :: AS_V

  ! unused..
  ! model_attenuation_variables
  !type model_attenuation_variables
  !  sequence
  !  double precision, dimension(N_SLS)            :: Qtau_s             ! tau_sigma
  !  integer :: dummy_pad ! padding 4 bytes to align the structure
  !  double precision, dimension(:), allocatable   :: Qmu                ! Shear Attenuation
  !  integer :: Qn                                                       ! Number of points
  !  double precision, dimension(:,:), allocatable :: Qtau_e             ! tau_epsilon
  !  double precision, dimension(:), allocatable   :: Qr                 ! Radius
  !  double precision, dimension(:), allocatable   :: QrDisc             ! Discontinuities Defined
  !  double precision, dimension(:), allocatable   :: Qomsb, Qomsb2      ! one_minus_sum_beta
  !  double precision, dimension(:,:), allocatable :: Qfc, Qfc2          ! factor_common
  !  double precision, dimension(:), allocatable   :: Qsf, Qsf2          ! scale_factor
  !  integer, dimension(:), allocatable            :: Qrmin              ! Max and Mins of idoubling
  !  integer, dimension(:), allocatable            :: Qrmax              ! Max and Mins of idoubling
  !  integer, dimension(:), allocatable            :: interval_Q         ! Steps
  !end type model_attenuation_variables
  !type (model_attenuation_variables) :: AM_V
  ! model_attenuation_variables

  end module model_attenuation_par

!
!------------------------------------------------------------------------
!

  subroutine model_attenuation_broadcast()

! standard routine to setup model

  use constants, only: N_SLS,myrank

  use shared_parameters, only: ATT_F_C_SOURCE
  use regions_mesh_par2, only: tau_s_store

  implicit none

  ! main process determines period ranges
  if (myrank == 0) call read_attenuation_model()

  ! broadcasts to all others
  call bcast_all_singledp(ATT_F_C_SOURCE)
  call bcast_all_dp(tau_s_store,N_SLS)

  end subroutine model_attenuation_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_attenuation_model()

  use constants, only: N_SLS

  use shared_parameters, only: MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,ATT_F_C_SOURCE

  use regions_mesh_par2, only: tau_s_store

  implicit none

  ! local parameters
  double precision :: f_c_source,min_period,max_period
  double precision,dimension(N_SLS) :: tau_s

  ! initializes
  tau_s(:) = 0.d0

  ! attenuation period range
  min_period = MIN_ATTENUATION_PERIOD
  max_period = MAX_ATTENUATION_PERIOD

  ! gets stress relaxation times
  call attenuation_tau_sigma(tau_s, N_SLS, min_period, max_period)

  ! gets center frequency (f_c_source) of simulation period band
  call attenuation_source_frequency(f_c_source, min_period, max_period)

  ! stores values
  tau_s_store(:) = tau_s(:)
  ATT_F_C_SOURCE = f_c_source

  end subroutine read_attenuation_model


!
!-------------------------------------------------------------------------------------------------
!

! This subroutine is hackish.  It could probably all be moved to an input attenuation file.
! Actually all the velocities, densities and attenuations could be moved to separate input
! files rather than be defined within the CODE
!
! All this subroutine does is define the Attenuation vs Radius and then Compute the Attenuation
! Variables (tau_sigma and tau_epsilon ( or tau_mu) )
  subroutine model_attenuation_setup(REFERENCE_1D_MODEL,CRUSTAL)

  use constants

  use regions_mesh_par2, only: tau_s_store

  use model_1dref_par, only: &
    NR_REF,Mref_V_Qmu_ref ! Mref_V_radius_ref

  use model_ak135_par, only: &
    NR_AK135F_NO_MUD,Mak135_V_Qmu_ak135 ! Mak135_V_radius_ak135

  use model_1066a_par, only: &
    NR_1066A,M1066a_V_Qmu_1066a ! M1066a_V_radius_1066a

  use model_sea1d_par, only: &
    NR_SEA1D,SEA1DM_V_Qmu_sea1d ! SEA1DM_V_radius_sea1d

  use model_case65tay_par, only: &
    NR_case65TAY,Mcase65TAY_V_Qmu

  use model_vpremoon_par, only: &
    NR_VPREMOON_layers,VPREMOON_Qmu_original

  implicit none

  integer,intent(in) :: REFERENCE_1D_MODEL
  logical,intent(in) :: CRUSTAL

  ! local parameters
  double precision :: tau_e(N_SLS)
  ! shear Attenuation
  double precision, dimension(:), allocatable   :: Qmu
  ! number of points
  integer :: Qn
  integer :: i,ier

  ! parameter definitions
  !double precision, parameter :: R120 = 6251.d3 ! as defined by IASP91

  ! initializes
  Qn = 0

  ! uses "pure" 1D models including their 1D-crust profiles
  ! (uses USE_EXTERNAL_CRUSTAL_MODEL set to false)
  select case(REFERENCE_1D_MODEL)
  case (REFERENCE_MODEL_PREM, &
        REFERENCE_MODEL_IASP91, &
        REFERENCE_MODEL_JP1D)
    ! PREM Q layers
    Qn = 12

  case (REFERENCE_MODEL_AK135F_NO_MUD)
    ! redefines "pure" 1D model without crustal modification
    call define_model_ak135(.false.)
    Qn = NR_AK135F_NO_MUD

  case (REFERENCE_MODEL_1066A)
    ! redefines "pure" 1D model without crustal modification
    call define_model_1066a(.false.)
    Qn = NR_1066A

  case (REFERENCE_MODEL_1DREF)
    ! redefines "pure" 1D model without crustal modification
    call define_model_1dref(.false.)
    Qn = NR_REF

  case (REFERENCE_MODEL_SEA1D)
    ! redefines "pure" 1D model without crustal modification
    call define_model_sea1d(.false.)
    Qn = NR_SEA1D

  case (REFERENCE_MODEL_SOHL, &
        REFERENCE_MODEL_SOHL_B)
    ! Mars, Q from PREM
    Qn = 12

  case (REFERENCE_MODEL_CASE65TAY)
    ! redefines "pure" 1D model without crustal modification
    call define_model_case65TAY(.false.)
    Qn = NR_case65TAY

  case (REFERENCE_MODEL_VPREMOON)
    ! Moon
    ! no need to redefine 1D model, Qmu_original array contains original values (without CRUSTAL modification)
    Qn = NR_VPREMOON_layers

  case default
    call exit_MPI(myrank, 'Error attenuation setup: Reference 1D Model Not recognized')
  end select

  ! checks
  if (Qn < 1) call exit_MPI(myrank,'Error Qmu array has size zero')

  ! sets up attenuation storage (for all possible Qmu values defined in the 1D models)
  allocate(Qmu(Qn),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary Qmu arrays')
  Qmu(:) = 0.d0

  select case(REFERENCE_1D_MODEL)
  case (REFERENCE_MODEL_PREM, &
        REFERENCE_MODEL_IASP91, &
        REFERENCE_MODEL_JP1D)
    ! PREM Q values
    ! radius = (/   0.0d0,    RICB,  RICB,  RCMB,    RCMB,    R670,    R670,   R220,    R220,    R80,     R80, R_PLANET /)
    Qmu(:) = (/  84.6d0,  84.6d0, 0.0d0, 0.0d0, 312.0d0, 312.0d0, 143.0d0, 143.0d0, 80.0d0, 80.0d0, 600.0d0, 600.0d0 /)

  case (REFERENCE_MODEL_AK135F_NO_MUD)
    ! radius = Mak135_V_radius_ak135(:)
    Qmu(:) = Mak135_V_Qmu_ak135(:)

  case (REFERENCE_MODEL_1066A)
    ! radius = M1066a_V_radius_1066a(:)
    Qmu(:) = M1066a_V_Qmu_1066a(:)

  case (REFERENCE_MODEL_1DREF)
    ! radius = Mref_V_radius_ref(:)
    Qmu(:) = Mref_V_Qmu_ref(:)

  case (REFERENCE_MODEL_SEA1D)
    ! radius = SEA1DM_V_radius_sea1d(:)
    Qmu(:) = SEA1DM_V_Qmu_sea1d(:)

  case (REFERENCE_MODEL_SOHL, &
        REFERENCE_MODEL_SOHL_B)
    ! Mars: todo - future use attenuation model by Nimmo & Faul, 2013?
    !       at the moment, this is taken from PREM
    ! radius = (/   0.0d0,    RICB,  RICB,  RCMB,    RCMB,    R670, R670,   R220,    R220,    R80,     R80, R_PLANET /)
    Qmu(:) = (/  84.6d0,  84.6d0, 0.0d0, 0.0d0, 312.0d0, 312.0d0, 143.0d0, 143.0d0, 80.0d0, 80.0d0, 600.0d0, 600.0d0 /)

  case (REFERENCE_MODEL_CASE65TAY)
    ! values as defined in model
    Qmu(:) = Mcase65TAY_V_Qmu(:)

  case (REFERENCE_MODEL_VPREMOON)
    ! Moon
    Qmu(:) = VPREMOON_Qmu_original(:)

  case default
    call exit_MPI(myrank, 'Error attenuation setup: Reference 1D Model values missing')
  end select

  ! sets up storage of relaxation times
  do i = 1,Qn
    call model_attenuation_getstored_tau(Qmu(i), tau_s_store, tau_e)
  enddo

  ! re-defines 1D models with crustal modification if necessary
  if (CRUSTAL) then
    select case(REFERENCE_1D_MODEL)
    case (REFERENCE_MODEL_AK135F_NO_MUD)
      ! redefines 1D model with crustal modification
      call define_model_ak135(CRUSTAL)

    case (REFERENCE_MODEL_1066A)
      ! redefines 1D model with crustal modification
      call define_model_1066a(CRUSTAL)

    case (REFERENCE_MODEL_1DREF)
      ! redefines 1D model with crustal modification
      call define_model_1dref(CRUSTAL)

    case (REFERENCE_MODEL_SEA1D)
      ! redefines 1D model with crustal modification
      call define_model_sea1d(CRUSTAL)

    case (REFERENCE_MODEL_CASE65TAY)
      ! redefines 1D model with crustal modification
      call define_model_case65TAY(CRUSTAL)
    end select
  endif

  end subroutine model_attenuation_setup

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_attenuation_getstored_tau(Qmu_in, tau_s, tau_e)

! includes min_period, max_period, and N_SLS

  use constants
  use shared_parameters, only: MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD

  implicit none

  double precision,intent(inout) :: Qmu_in
  double precision, dimension(N_SLS),intent(inout) :: tau_s, tau_e

  ! local parameters
  integer :: rw
  double precision :: min_period,max_period

  ! attenuation period band
  min_period = MIN_ATTENUATION_PERIOD
  max_period = MAX_ATTENUATION_PERIOD

  ! tries to retrieve value from storage table
  ! READ
  rw = 1
  call model_attenuation_storage(Qmu_in, tau_e, rw)
  ! checks if read was successful
  if (rw > 0) return

  ! value not found in table, will need to create appropriate entries
  call attenuation_invert_by_simplex(min_period, max_period, N_SLS, Qmu_in, tau_s, tau_e)

  ! stores value in storage table
  ! WRITE
  rw = -1
  call model_attenuation_storage(Qmu_in, tau_e, rw)

  end subroutine model_attenuation_getstored_tau

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_attenuation_storage(Qmu, tau_e, rw)

  use constants

  use model_attenuation_par, only: AM_S

  implicit none

  double precision,intent(inout) :: Qmu
  double precision, dimension(N_SLS),intent(inout) :: tau_e
  integer,intent(inout) :: rw

  ! local parameters
  double precision :: Qmu_new
  integer :: num_total_Q,Qtmp,ier
  integer, save :: first_time_called = 1

  if (first_time_called == 1) then
    first_time_called = 0
    AM_S%Q_resolution = 10**ATTENUATION_COMP_RESOLUTION

    ! total number of Q values for storage table
    num_total_Q = AM_S%Q_resolution * ATTENUATION_COMP_MAXIMUM

    ! allocates Q table for precomputing relaxation times
    allocate(AM_S%tau_e_storage(N_SLS, num_total_Q), &
             AM_S%Qmu_storage(num_total_Q),stat=ier)
    if (ier /= 0) stop 'Error allocating AM_S arrays'
    AM_S%tau_e_storage(:,:) = 0.d0
    AM_S%Qmu_storage(:) = -1
  endif

  ! limits shear attenuation value
  ! by default, value set in constants.h for ATTENUATION_COMP_MAXIMUM = 9000
  ! (higher Q -> attenuation becomes negligible)
  !
  ! zero/negative Q values or larger Q values will be limited to this maximum values
  if (Qmu <= 0.0d0) Qmu = ATTENUATION_COMP_MAXIMUM
  if (Qmu > ATTENUATION_COMP_MAXIMUM) Qmu = ATTENUATION_COMP_MAXIMUM

  ! won't get executed since Qmu not zero anymore, but left here...
  if (rw > 0 .and. Qmu == 0.0d0) then
    ! read
    Qmu = 0.0d0
    tau_e(:) = 0.0d0
    return
  endif

  ! Generate index for Storage Array
  ! and Recast Qmu using this index

  ! by default: resolution is Q_resolution = 10
  ! converts Qmu to an array integer index:
  ! e.g. Qmu = 150.31 -> Qtmp = 150.31 * 10 = int( 1503.10 ) = 1503
  Qtmp    = int(Qmu * dble(AM_S%Q_resolution))

  ! limits
  if (Qtmp < 1) Qtmp = 1
  if (Qtmp > AM_S%Q_resolution * ATTENUATION_COMP_MAXIMUM) Qtmp = AM_S%Q_resolution * ATTENUATION_COMP_MAXIMUM

  ! rounds to corresponding double value:
  ! e.g. Qmu_new = dble( 1503 ) / dble(10) = 150.30
  ! but Qmu_new is not used any further...
  Qmu_new = dble(Qtmp) / dble(AM_S%Q_resolution)

  if (rw > 0) then
    ! READ
    if (AM_S%Qmu_storage(Qtmp) > 0) then
      ! READ SUCCESSFUL
      tau_e(:)   = AM_S%tau_e_storage(:, Qtmp)
      Qmu        = AM_S%Qmu_storage(Qtmp)
      rw = 1
    else
      ! READ NOT SUCCESSFUL
      rw = -1
    endif
  else
    ! WRITE SUCCESSFUL
    AM_S%tau_e_storage(:,Qtmp)    = tau_e(:)
    AM_S%Qmu_storage(Qtmp)        = Qmu
    rw = 1
  endif

  end subroutine model_attenuation_storage


!
!-------------------------------------------------------------------------------------------------
!

  subroutine attenuation_source_frequency(f_c_source, min_period, max_period)

! Determine the Source Frequency

  use constants, only: myrank

  implicit none

  double precision,intent(out) :: f_c_source
  double precision,intent(in) :: min_period, max_period

  ! local parameters
  double precision :: f1, f2
  double precision, parameter :: TOL_ZERO = 1.d-30

  ! min and max frequencies of the simulation
  f1 = 1.0d0 / max_period
  f2 = 1.0d0 / min_period

  ! takes the frequency at the logarithmic mean between min and max
  ! note: previous versions had a scaling factor 1000.d0 added and removed again when read in the solver.
  f_c_source = 10.0d0**(0.5 * (log10(f1) + log10(f2)))

  ! checks
  if (f_c_source < TOL_ZERO .or. f_c_source /= f_c_source) then
    print *,'Error: invalid center frequency ',f_c_source,'in attenuation_source_frequency() routine'
    call exit_MPI(myrank,'Error invalid center frequency')
  endif

  end subroutine attenuation_source_frequency

!
!-------------------------------------------------------------------------------------------------
!

  subroutine attenuation_tau_sigma(tau_s, n, min_period, max_period)

! Set the Tau_sigma (tau_s) to be equally spaced in log10 frequency

  use constants, only: PI,myrank

  implicit none

  integer,intent(in) :: n
  double precision,intent(inout) :: tau_s(n)
  double precision,intent(in) :: min_period, max_period

  ! local parameters
  double precision :: f1, f2
  double precision :: exp1, exp2
  double precision :: dexpval,fac
  double precision, parameter :: TOL_ZERO = 1.d-30
  integer :: i

  ! checks
  if (max_period < TOL_ZERO) &
    call exit_MPI(myrank,'Error invalid maximum period in attenuation_tau_sigma(), cannot be zero or negative.')
  if (min_period < TOL_ZERO) &
    call exit_MPI(myrank,'Error invalid maximum period in attenuation_tau_sigma(), cannot be zero or negative.')

  ! min and max frequencies
  f1 = 1.0d0 / max_period
  f2 = 1.0d0 / min_period

  ! log(f)
  exp1 = log10(f1)
  exp2 = log10(f2)

  ! increment in log-space
  dexpval = (exp2-exp1) / ((n*1.0d0) - 1)

  ! equally distributed in log-space
  do i = 1,n
    fac = PI * 2.0d0 * 10**(exp1 + (i - 1) * 1.0d0 * dexpval)
    ! checks
    if (abs(fac) <= TOL_ZERO) then
      print *,'Error: invalid denominator ',fac,' for sls ',i,' in attenuation_tau_sigma() routine'
      stop 'Error invalid denominator in attenuation_tau_sigma() routine'
    endif
    ! relaxation time
    tau_s(i) = 1.d0 / fac
  enddo

  end subroutine attenuation_tau_sigma

!
!-------------------------------------------------------------------------------------------------
!

  subroutine attenuation_invert_by_simplex(t2, t1, n, Q_real, tau_s, tau_e)

  use constants, only: PI,myrank

  implicit none

  ! Input / Output
  double precision,intent(in) :: t1, t2  ! min/max of attenuation period band
  double precision,intent(in) :: Q_real
  integer,intent(in) :: n
  double precision, dimension(n),intent(inout) :: tau_s, tau_e

  ! local parameters
  integer :: i, iterations, err,prnt
  double precision :: f1, f2, exp1,exp2,dexpval, min_value

  ! frequencies (to compare calculated and target results)
  integer, parameter :: nf = 100
  double precision, dimension(nf) :: f

  double precision, external :: attenuation_eval

  ! Values to be passed into the simplex minimization routine
  iterations = -1
  min_value  = -1.0e-4
  err        = 0
  prnt       = 0

  ! Determine the min and max frequencies
  f1 = 1.0d0 / t1
  f2 = 1.0d0 / t2

  ! Determine the exponents of the frequencies
  exp1 = log10(f1)
  exp2 = log10(f2)

  if (f2 < f1 .or. Q_real < 0.0d0 .or. n < 1) then
    call exit_MPI(myrank, 'frequencies flipped or Q less than zero or N_SLS < 0')
  endif

  ! Determine the Frequencies at which to compare solutions
  ! The frequencies should be equally spaced in log10 frequency
  do i = 1,nf
     f(i) = exp1 + ((i-1)*1.0d0 * (exp2-exp1) / ((nf-1)*1.0d0))
  enddo

  ! Set the Tau_sigma (tau_s) to be equally spaced in log10 frequency
  dexpval = (exp2-exp1) / ((n*1.0d0) - 1)
  do i = 1,n
     tau_s(i) = 1.0 / (PI * 2.0d0 * 10**(exp1 + (i - 1)* 1.0d0 *dexpval))
  enddo

  ! note: for frequency dependent Q models, one would need to define Q(f) at the frequencies f(i) defined above.
  !       this would require to setup the arrays AS_V%Q and AS_V%iQ to have Q values for each frequency f_i.
  !
  !       here we assume a frequency-independent Q model, thus we only need to specify single values
  !       in AS_V%Q and AS_V%iQ, respectively, based Q_real.
  !
  !       this will be done in the following routine call.
  !
  ! Shove the parameters into the module
  call attenuation_simplex_setup(nf,n,f,Q_real,tau_s)

  ! Set the Tau_epsilon (tau_e) to an initial value at omega*tau = 1
  ! tan_delta = 1/Q = (tau_e - tau_s)/(2 * sqrt(tau e*tau_s))
  !    if we assume tau_e =~ tau_s
  !    we get the equation below
  do i = 1,n
     tau_e(i) = tau_s(i) + (tau_s(i) * 2.0d0/Q_real)
  enddo

  ! Run a simplex search to determine the optimum values of tau_e
  call fminsearch(attenuation_eval, tau_e, n, iterations, min_value, prnt, err)

  if (err > 0) then
     write(*,*)'Search did not converge for an attenuation of ', Q_real
     write(*,*)'    Iterations: ', iterations
     write(*,*)'    Min Value:  ', min_value
     write(*,*)'    Aborting program'
     call exit_MPI(myrank,'attenuation_simplex: Search for Strain relaxation times did not converge')
  endif

  call attenuation_simplex_finish()

  end subroutine attenuation_invert_by_simplex

!
!-------------------------------------------------------------------------------------------------
!

  subroutine attenuation_simplex_finish()

  use model_attenuation_par, only: AS_V

  implicit none

  deallocate(AS_V%f)
  deallocate(AS_V%tau_s)

  end subroutine attenuation_simplex_finish


!
!-------------------------------------------------------------------------------------------------
!

!   - Inserts necessary parameters into the module attenuation_simplex_variables
!   - See module for explanation
  subroutine attenuation_simplex_setup(nf_in, nsls_in, f_in, Q_in, tau_s_in)

  use model_attenuation_par, only: AS_V

  implicit none

  integer,intent(in) :: nf_in, nsls_in
  double precision,intent(in) :: Q_in
  double precision, dimension(nf_in),intent(in) :: f_in
  double precision, dimension(nsls_in),intent(in) :: tau_s_in

  ! local parameters
  integer :: ier

  ! allocates structure arrays
  allocate(AS_V%f(nf_in), &
           AS_V%tau_s(nsls_in),stat=ier)
  if (ier /= 0) stop 'Error allocating AS_V arrays'

  ! initializes
  AS_V%nf       = nf_in
  AS_V%nsls     = nsls_in
  AS_V%f(:)     = f_in(:)
  AS_V%Q        = Q_in
  AS_V%iQ       = 1.0d0/AS_V%Q
  AS_V%tau_s(:) = tau_s_in(:)

  end subroutine attenuation_simplex_setup

!
!-------------------------------------------------------------------------------------------------
!

!   - Computes the Moduli (Maxwell Solid) for a series of
!         Standard Linear Solids
!   - Computes M1 and M2 parameters after Dahlen and Tromp pp.203
!         here called B and A after Liu et al. 1976
!   - Another formulation uses Kelvin-Voigt Solids and computes
!         Compliences J1 and J2 after Dahlen and Tromp pp.203
!
!   Input
!     nf    = Number of Frequencies
!     nsls  = Number of Standard Linear Solids
!     f     = Frequencies (in log10 of frequencies)
!                dimension(nf)
!     tau_s = Tau_sigma  Stress relaxation time (see References)
!                dimension(nsls)
!     tau_e = Tau_epsilon Strain relaxation time (see References)
!                dimension(nsls)!
!   Output
!     B     = Real Moduli      ( M2 Dahlen and Tromp pp.203 )
!                dimension(nf)
!     A     = Imaginary Moduli ( M1 Dahlen and Tromp pp.203 )
!                dimension(nf)
!
!   Dahlen and Tromp, 1998
!      Theoretical Global Seismology
!
!   Liu et al. 1976
!      Velocity dispersion due to anelasticity: implications for seismology and mantle composition
!      Geophys, J. R. asts. Soc, Vol 47, pp. 41-58
  subroutine attenuation_maxwell(nf,nsls,f,tau_s,tau_e,B,A)

  use constants

  implicit none

  ! Input
  integer,intent(in) :: nf, nsls
  double precision, dimension(nf),intent(in)   :: f
  double precision, dimension(nsls),intent(in) :: tau_s, tau_e
  ! Output
  double precision, dimension(nf),intent(out)   :: A,B

  integer :: i,j
  double precision :: w, demon

  A(:) = 1.0d0 -  nsls*1.0d0
  B(:) = 0.0d0
  do i = 1,nf
    ! angular frequency omega
    w = 2.0d0 * PI * 10**f(i)
    do j = 1,nsls
      ! denominator
      demon = 1.0d0 + w**2 * tau_s(j)**2
      ! moduli
      A(i) = A(i) + ((1.0d0 + (w**2 * tau_e(j) * tau_s(j)))/ demon)
      B(i) = B(i) + ((w * (tau_e(j) - tau_s(j))) / demon)
    enddo
  enddo

  end subroutine attenuation_maxwell

!
!-------------------------------------------------------------------------------------------------
!

!    - Computes the misfit from a set of relaxation parameters
!          given a set of frequencies and target attenuation
!    - Evaluates only at the given frequencies
!    - Evaluation is done with an L2 norm
!
!    Input
!      Xin = Tau_epsilon, Strain Relaxation Time
!                Note: Tau_sigma the Stress Relaxation Time is loaded
!                      with attenuation_simplex_setup and stored in
!                      attenuation_simplex_variables
!
!    Xi = Sum_i^N sqrt [ (1/Qc_i - 1/Qt_i)^2 / 1/Qt_i^2 ]
!
!     where Qc_i is the computed attenuation at a specific frequency
!           Qt_i is the desired attenuation at that frequency
!
!    Uses attenuation_simplex_variables to store constant values
!
!    See attenuation_simplex_setup
!
  double precision function attenuation_eval(Xin)

  use constants, only: HUGEVAL

  use model_attenuation_par, only: AS_V

  implicit none

  ! Input
  double precision, dimension(AS_V%nsls),intent(in) :: Xin

  ! local parameters
  double precision, dimension(AS_V%nsls) :: tau_e
  double precision, dimension(AS_V%nf)   :: A, B, tan_delta
  integer :: i
  double precision :: xi, iQ2

  tau_e = Xin

  call attenuation_maxwell(AS_V%nf,AS_V%nsls,AS_V%f,AS_V%tau_s,tau_e,B,A)

  do i = 1,AS_V%nf
    if (abs(A(i)) > 0.0d0) then
      tan_delta(i) = B(i) / A(i)
    else
      tan_delta(i) = HUGEVAL
    endif
  enddo

  attenuation_eval = 0.0d0
  iQ2 = AS_V%iQ**2
  do i = 1,AS_V%nf
     xi = sqrt(( ( (tan_delta(i) - AS_V%iQ) ** 2 ) / iQ2 ))
     attenuation_eval = attenuation_eval + xi
  enddo

  end function attenuation_eval

!
!-------------------------------------------------------------------------------------------------
!

! subroutine fminsearch
!   - Computes the minimization of funk(x(n)) using the simplex method
!   - This subroutine is copied from Matlab fminsearch.m
!         and modified to suit my nefarious needs
!   Input
!     funk = double precision function with one input parameter
!                double precision function the_funk(x)
!     x    = Input/Output
!               variables to be minimized
!               dimension(n)
!            Input:  Initial Value
!            Output: Minimized Value
!     n    = number of variables
!     itercount = Input/Output
!                 Input:  maximum number of iterations
!                         if < 0 default is used (200 * n)
!                 Output: total number of iterations on output
!     tolf      = Input/Output
!                 Input:  minimum tolerance of the function funk(x)
!                 Output: minimum value of funk(x)(i.e. "a" solution)
!     prnt      = Input
!                 3 => report every iteration
!                 4 => report every iteration, total simplex
!     err       = Output
!                 0 => Normal execution, converged within desired range
!                 1 => function evaluation exceeded limit
!                 2 => Iterations exceeded limit
!
!     See Matlab fminsearch
  subroutine fminsearch(funk, x, n, itercount, tolf, prnt, err)

  implicit none

  ! Input
  double precision, external :: funk

  integer,intent(in) :: n
  double precision,intent(inout) :: x(n) ! Also Output
  integer,intent(inout) :: itercount, err
  integer,intent(in) :: prnt
  double precision,intent(inout) :: tolf

  !Internal
  integer :: i,j, how
  integer, parameter :: none             = 0
  integer, parameter :: initial          = 1
  integer, parameter :: expand           = 2
  integer, parameter :: reflect          = 3
  integer, parameter :: contract_outside = 4
  integer, parameter :: contract_inside  = 5
  integer, parameter :: shrink           = 6

  integer :: maxiter, maxfun
  integer :: func_evals
  double precision :: tolx

  double precision :: rho, chi, psi, sigma
  double precision :: xin(n), y(n), v(n,n+1), fv(n+1)
  double precision :: vtmp(n,n+1)
  double precision :: usual_delta, zero_term_delta
  double precision :: xbar(n), xr(n), fxr, xe(n), fxe, xc(n), fxc, fxcc, xcc(n)
  integer :: place(n+1)

  double precision,external :: max_size_simplex, max_value

  rho   = 1.0d0
  chi   = 2.0d0
  psi   = 0.5d0
  sigma = 0.5d0

  if (itercount > 0) then
     maxiter = itercount
  else
     maxiter = 200 * n
  endif
  itercount = 0
  maxfun  = 200 * n

  if (tolf > 0.0d0) then
     tolx = 1.0e-4
  else
     tolx = 1.0e-4
     tolf = 1.0e-4
  endif

  err = 0

  xin    = x
  v(:,:) = 0.0d0
  fv(:)  = 0.0d0

  v(:,1) = xin
  x      = xin

  fv(1) = funk(xin)

  usual_delta = 0.05
  zero_term_delta = 0.00025

  do j = 1,n
     y = xin
     if (y(j) /= 0.0d0) then
        y(j) = (1.0d0 + usual_delta) * y(j)
     else
        y(j) = zero_term_delta
     endif
     v(:,j+1) = y
     x(:) = y
     fv(j+1) = funk(x)
  enddo

  call heap_sort_local(n+1, fv, place)

  do i = 1,n+1
    vtmp(:,i) = v(:,place(i))
  enddo
  v = vtmp

  how = initial
  itercount = 1
  func_evals = n+1
  if (prnt == 3) then
     write(*,*)'Iterations   Funk Evals   Value How'
     write(*,*)itercount, func_evals, fv(1), how
  endif
  if (prnt == 4) then
     write(*,*)'How: ',how
     write(*,*)'V: ', v
     write(*,*)'fv: ',fv
     write(*,*)'evals: ',func_evals
  endif

  do while (func_evals < maxfun .and. itercount < maxiter)

     if (max_size_simplex(v,n) <= tolx .and. max_value(fv,n+1) <= tolf) then
       goto 666
     endif
     how = none

     ! xbar = average of the n (NOT n+1) best points
     !     xbar = sum(v(:,1:n), 2)/n
     xbar(:) = 0.0d0
     do i = 1,n
        do j = 1,n
           xbar(i) = xbar(i) + v(i,j)
        enddo
        xbar(i) = xbar(i) / (n*1.0d0)
     enddo
     xr = (1 + rho)*xbar - rho*v(:,n+1)
     x(:) = xr
     fxr = funk(x)
     func_evals = func_evals + 1
     if (fxr < fv(1)) then
        ! Calculate the expansion point
        xe = (1 + rho*chi)*xbar - rho*chi*v(:,n+1)
        x = xe
        fxe = funk(x)
        func_evals = func_evals+1
        if (fxe < fxr) then
           v(:,n+1) = xe
           fv(n+1) = fxe
           how = expand
        else
           v(:,n+1) = xr
           fv(n+1) = fxr
           how = reflect
        endif
     else ! fv(:,1) <= fxr
        if (fxr < fv(n)) then
           v(:,n+1) = xr
           fv(n+1) = fxr
           how = reflect
        else ! fxr >= fv(:,n)
           ! Perform contraction
           if (fxr < fv(n+1)) then
              ! Perform an outside contraction
              xc = (1 + psi*rho)*xbar - psi*rho*v(:,n+1)
              x(:) = xc
              fxc = funk(x)
              func_evals = func_evals+1

              if (fxc <= fxr) then
                 v(:,n+1) = xc
                 fv(n+1) = fxc
                 how = contract_outside
              else
                 ! perform a shrink
                 how = shrink
              endif
           else
              ! Perform an inside contraction
              xcc = (1-psi)*xbar + psi*v(:,n+1)
              x(:) = xcc
              fxcc = funk(x)
              func_evals = func_evals+1

              if (fxcc < fv(n+1)) then
                 v(:,n+1) = xcc
                 fv(n+1) = fxcc
                 how = contract_inside
              else
                 ! perform a shrink
                 how = shrink
              endif
           endif
           if (how == shrink) then
              do j=2,n+1
                 v(:,j)=v(:,1)+sigma*(v(:,j) - v(:,1))
                 x(:) = v(:,j)
                 fv(j) = funk(x)
              enddo
              func_evals = func_evals + n
           endif
        endif
     endif

     call heap_sort_local(n+1, fv, place)

     do i = 1,n+1
       vtmp(:,i) = v(:,place(i))
     enddo
     v = vtmp

     itercount = itercount + 1
     if (prnt == 3) then
        write(*,*)itercount, func_evals, fv(1), how
     else if (prnt == 4) then
        write(*,*)
        write(*,*)'How: ',how
        write(*,*)'v: ',v
        write(*,*)'fv: ',fv
        write(*,*)'evals: ',func_evals
     endif
  enddo

  if (func_evals > maxfun) then
     write(*,*)'function evaluations exceeded prescribed limit', maxfun
     err = 1
  endif
  if (itercount > maxiter) then
     write(*,*)'iterations exceeded prescribed limit', maxiter
     err = 2
  endif

666 continue
  x = v(:,1)
  tolf = fv(1)

  end subroutine fminsearch

!
!-------------------------------------------------------------------------------------------------
!

!    - Finds the maximum value of the difference of between the first
!          value and the remaining values of a vector
!    Input
!      fv = Input
!             Vector
!             dimension(n)
!      n  = Input
!             Length of fv
!
!      returns:
!         Xi = max( || fv(1)- fv(i) || ) for i=2:n
!
  double precision function max_value(fv,n)
  implicit none
  integer,intent(in) :: n
  double precision, intent(in) :: fv(n)

  integer :: i
  double precision :: m, z

  m = 0.0d0
  do i = 2,n
     z = abs(fv(1) - fv(i))
     if (z > m) then
        m = z
     endif
  enddo

  max_value = m

  end function max_value

!
!-------------------------------------------------------------------------------------------------
!

!   - Determines the maximum distance between two point in a simplex
!   Input
!     v  = Input
!            Simplex Vertices
!            dimension(n, n+1)
!     n  = Pseudo Length of n
!
!     returns:
!       Xi = max( max( || v(:,1) - v(:,i) || ) ) for i=2:n+1
!
  double precision function max_size_simplex(v,n)
  implicit none
  integer, intent(in) :: n
  double precision, intent(in) :: v(n,n+1)

  integer :: i,j
  double precision :: m, z

  m = 0.0d0
  do i = 1,n
     do j = 2,n+1
        z = abs(v(i,j) - v(i,1))
        if (z > m) then
           m = z
        endif
     enddo
  enddo

  max_size_simplex = m

  end function max_size_simplex

