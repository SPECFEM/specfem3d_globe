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
!
!  This portion of the SPECFEM3D Code was written by:
!  Brian Savage while at
!     California Institute of Technology
!     Department of Terrestrial Magnetism / Carnegie Institute of Washington
!     University of Rhode Island
!
!  <savage@uri.edu>.
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

  subroutine model_attenuation_broadcast(myrank,AM_V,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD)

! standard routine to setup model

  use constants

  implicit none

! model_attenuation_variables
  type model_attenuation_variables
    sequence
    double precision min_period, max_period
    double precision                          :: QT_c_source        ! Source Frequency
    double precision, dimension(:), pointer   :: Qtau_s             ! tau_sigma
    double precision, dimension(:), pointer   :: QrDisc             ! Discontinuities Defined
    double precision, dimension(:), pointer   :: Qr                 ! Radius
    double precision, dimension(:), pointer   :: Qmu                ! Shear Attenuation
    double precision, dimension(:,:), pointer :: Qtau_e             ! tau_epsilon
    double precision, dimension(:), pointer   :: Qomsb, Qomsb2      ! one_minus_sum_beta
    double precision, dimension(:,:), pointer :: Qfc, Qfc2          ! factor_common
    double precision, dimension(:), pointer   :: Qsf, Qsf2          ! scale_factor
    integer, dimension(:), pointer            :: Qrmin              ! Max and Mins of idoubling
    integer, dimension(:), pointer            :: Qrmax              ! Max and Mins of idoubling
    integer, dimension(:), pointer            :: interval_Q                 ! Steps
    integer                                   :: Qn                 ! Number of points
    integer dummy_pad ! padding 4 bytes to align the structure
  end type model_attenuation_variables

  type (model_attenuation_variables) AM_V
! model_attenuation_variables

  integer :: MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD
  integer :: myrank
  integer :: ier

  allocate(AM_V%Qtau_s(N_SLS),stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating Qtau_s array')

  ! master process determines period ranges
  if (myrank == 0) call read_attenuation_model(MIN_ATTENUATION_PERIOD, MAX_ATTENUATION_PERIOD, AM_V)

  ! broadcasts to all others
  call bcast_all_singledp(AM_V%min_period)
  call bcast_all_singledp(AM_V%max_period)
  call bcast_all_singledp(AM_V%QT_c_source)
  call bcast_all_dp(AM_V%Qtau_s,   N_SLS)

  end subroutine model_attenuation_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_attenuation_model(min_att_period, max_att_period, AM_V)

  use constants

  implicit none

! model_attenuation_variables
  type model_attenuation_variables
    sequence
    double precision min_period, max_period
    double precision                          :: QT_c_source        ! Source Frequency
    double precision, dimension(:), pointer   :: Qtau_s             ! tau_sigma
    double precision, dimension(:), pointer   :: QrDisc             ! Discontinuities Defined
    double precision, dimension(:), pointer   :: Qr                 ! Radius
    double precision, dimension(:), pointer   :: Qmu                ! Shear Attenuation
    double precision, dimension(:,:), pointer :: Qtau_e             ! tau_epsilon
    double precision, dimension(:), pointer   :: Qomsb, Qomsb2      ! one_minus_sum_beta
    double precision, dimension(:,:), pointer :: Qfc, Qfc2          ! factor_common
    double precision, dimension(:), pointer   :: Qsf, Qsf2          ! scale_factor
    integer, dimension(:), pointer            :: Qrmin              ! Max and Mins of idoubling
    integer, dimension(:), pointer            :: Qrmax              ! Max and Mins of idoubling
    integer, dimension(:), pointer            :: interval_Q                 ! Steps
    integer                                   :: Qn                 ! Number of points
    integer dummy_pad ! padding 4 bytes to align the structure
  end type model_attenuation_variables

  type (model_attenuation_variables) AM_V
! model_attenuation_variables

  integer :: min_att_period, max_att_period

  AM_V%min_period = min_att_period * 1.0d0
  AM_V%max_period = max_att_period * 1.0d0

  call attenuation_tau_sigma(AM_V%Qtau_s, N_SLS, AM_V%min_period, AM_V%max_period)
  call attenuation_source_frequency(AM_V%QT_c_source, AM_V%min_period, AM_V%max_period)

  end subroutine read_attenuation_model


!
!-------------------------------------------------------------------------------------------------
!

! This Subroutine is Hackish.  It could probably all be moved to an input attenuation file.
! Actually all the velocities, densities and attenuations could be moved to separate input
! files rather than be defined within the CODE
!
! All this subroutine does is define the Attenuation vs Radius and then Compute the Attenuation
! Variables (tau_sigma and tau_epsilon ( or tau_mu) )
  subroutine model_attenuation_setup(myrank,REFERENCE_1D_MODEL,RICB,RCMB, &
                                    R670,R220,R80,AM_V,AM_S,AS_V,CRUSTAL)

  use constants

  use model_1dref_par, only: &
    NR_REF,Mref_V_radius_ref,Mref_V_Qmu_ref

  use model_ak135_par, only: &
    NR_AK135F_NO_MUD,Mak135_V_radius_ak135,Mak135_V_Qmu_ak135

  use model_1066a_par, only: &
    NR_1066A,M1066a_V_radius_1066a,M1066a_V_Qmu_1066a

  use model_sea1d_par, only: &
    NR_SEA1D,SEA1DM_V_radius_sea1d,SEA1DM_V_Qmu_sea1d

  implicit none

! model_attenuation_variables
  type model_attenuation_variables
    sequence
    double precision min_period, max_period
    double precision                          :: QT_c_source        ! Source Frequency
    double precision, dimension(:), pointer   :: Qtau_s             ! tau_sigma
    double precision, dimension(:), pointer   :: QrDisc             ! Discontinuities Defined
    double precision, dimension(:), pointer   :: Qr                 ! Radius
    double precision, dimension(:), pointer   :: Qmu                ! Shear Attenuation
    double precision, dimension(:,:), pointer :: Qtau_e             ! tau_epsilon
    double precision, dimension(:), pointer   :: Qomsb, Qomsb2      ! one_minus_sum_beta
    double precision, dimension(:,:), pointer :: Qfc, Qfc2          ! factor_common
    double precision, dimension(:), pointer   :: Qsf, Qsf2          ! scale_factor
    integer, dimension(:), pointer            :: Qrmin              ! Max and Mins of idoubling
    integer, dimension(:), pointer            :: Qrmax              ! Max and Mins of idoubling
    integer, dimension(:), pointer            :: interval_Q                 ! Steps
    integer                                   :: Qn                 ! Number of points
    integer dummy_pad ! padding 4 bytes to align the structure
  end type model_attenuation_variables

  type (model_attenuation_variables) AM_V
! model_attenuation_variables

! model_attenuation_storage_var
  type model_attenuation_storage_var
    sequence
    double precision, dimension(:,:), pointer :: tau_e_storage
    double precision, dimension(:), pointer :: Qmu_storage
    integer Q_resolution
    integer Q_max
  end type model_attenuation_storage_var

  type (model_attenuation_storage_var) AM_S
! model_attenuation_storage_var

! attenuation_simplex_variables
  type attenuation_simplex_variables
    sequence
    double precision Q  ! Q     = Desired Value of Attenuation or Q
    double precision iQ ! iQ    = 1/Q
    double precision, dimension(:), pointer ::  f
    ! f = Frequencies at which to evaluate the solution
    double precision, dimension(:), pointer :: tau_s
    ! tau_s = Tau_sigma defined by the frequency range and
    !             number of standard linear solids
    integer nf          ! nf    = Number of Frequencies
    integer nsls        ! nsls  = Number of Standard Linear Solids
  end type attenuation_simplex_variables

  type(attenuation_simplex_variables) AS_V
! attenuation_simplex_variables

  integer :: myrank,REFERENCE_1D_MODEL
  double precision :: RICB, RCMB, R670, R220, R80
  logical :: CRUSTAL

  ! local parameters
  double precision :: tau_e(N_SLS)
  integer :: i,ier

  ! parameter definitions
  double precision, parameter :: R120 = 6251.d3 ! as defined by IASP91

  ! uses "pure" 1D models including their 1D-crust profiles
  ! (uses USE_EXTERNAL_CRUSTAL_MODEL set to false)
  if (REFERENCE_1D_MODEL == REFERENCE_MODEL_PREM) then
    AM_V%Qn = 12
  else if (REFERENCE_1D_MODEL == REFERENCE_MODEL_IASP91) then
    AM_V%Qn = 12
  else if (REFERENCE_1D_MODEL == REFERENCE_MODEL_AK135F_NO_MUD) then
    ! redefines "pure" 1D model without crustal modification
    call define_model_ak135(.FALSE.)
    AM_V%Qn = NR_AK135F_NO_MUD
  else if (REFERENCE_1D_MODEL == REFERENCE_MODEL_1066A) then
    ! redefines "pure" 1D model without crustal modification
    call define_model_1066a(.FALSE.)
    AM_V%Qn = NR_1066A
  else if (REFERENCE_1D_MODEL == REFERENCE_MODEL_1DREF) then
    ! redefines "pure" 1D model without crustal modification
    call define_model_1dref(.FALSE.)
    AM_V%Qn = NR_REF
  else if (REFERENCE_1D_MODEL == REFERENCE_MODEL_JP1D) then
    AM_V%Qn = 12
  else if (REFERENCE_1D_MODEL == REFERENCE_MODEL_SEA1D) then
    ! redefines "pure" 1D model without crustal modification
    call define_model_sea1d(.FALSE.)
    AM_V%Qn = NR_SEA1D
  else
    call exit_MPI(myrank, 'Reference 1D Model Not recognized')
  endif

  ! sets up attenuation storage (for all possible Qmu values defined in the 1D models)
  allocate(AM_V%Qr(AM_V%Qn), &
           AM_V%Qmu(AM_V%Qn), &
           AM_V%interval_Q(AM_V%Qn), &
           AM_V%Qtau_e(N_SLS,AM_V%Qn), &
           stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating AM_V arrays')

  if (REFERENCE_1D_MODEL == REFERENCE_MODEL_PREM) then
    AM_V%Qr(:)     = (/    0.0d0,     RICB,  RICB,  RCMB,    RCMB,    R670,    R670,   R220,    R220,    R80,     R80, R_EARTH /)
    AM_V%Qmu(:)    = (/   84.6d0,   84.6d0, 0.0d0, 0.0d0, 312.0d0, 312.0d0, 143.0d0, 143.0d0, 80.0d0, 80.0d0, 600.0d0, 600.0d0 /)
  else if (REFERENCE_1D_MODEL == REFERENCE_MODEL_IASP91) then
    AM_V%Qr(:)     = (/    0.0d0,     RICB,  RICB,  RCMB,    RCMB,    R670,    R670,    R220,   R220,   R120,    R120, R_EARTH /)
    AM_V%Qmu(:)    = (/   84.6d0,   84.6d0, 0.0d0, 0.0d0, 312.0d0, 312.0d0, 143.0d0, 143.0d0, 80.0d0, 80.0d0, 600.0d0, 600.0d0 /)
  else if (REFERENCE_1D_MODEL == REFERENCE_MODEL_AK135F_NO_MUD) then
    AM_V%Qr(:)     = Mak135_V_radius_ak135(:)
    AM_V%Qmu(:)    = Mak135_V_Qmu_ak135(:)
  else if (REFERENCE_1D_MODEL == REFERENCE_MODEL_1066A) then
    AM_V%Qr(:)     = M1066a_V_radius_1066a(:)
    AM_V%Qmu(:)    = M1066a_V_Qmu_1066a(:)
  else if (REFERENCE_1D_MODEL == REFERENCE_MODEL_1DREF) then
    AM_V%Qr(:)     = Mref_V_radius_ref(:)
    AM_V%Qmu(:)    = Mref_V_Qmu_ref(:)
  else if (REFERENCE_1D_MODEL == REFERENCE_MODEL_JP1D) then
    AM_V%Qr(:)     = (/    0.0d0,     RICB,  RICB,  RCMB,    RCMB,    R670,    R670,    R220,   R220,   R120,    R120, R_EARTH /)
    AM_V%Qmu(:)    = (/   84.6d0,   84.6d0, 0.0d0, 0.0d0, 312.0d0, 312.0d0, 143.0d0, 143.0d0, 80.0d0, 80.0d0, 600.0d0, 600.0d0 /)
  else if (REFERENCE_1D_MODEL == REFERENCE_MODEL_SEA1D) then
    AM_V%Qr(:)     = SEA1DM_V_radius_sea1d(:)
    AM_V%Qmu(:)    = SEA1DM_V_Qmu_sea1d(:)
  endif

  do i = 1, AM_V%Qn
    call model_attenuation_getstored_tau(AM_V%Qmu(i), AM_V%QT_c_source, AM_V%Qtau_s, tau_e, AM_V, AM_S,AS_V)
    AM_V%Qtau_e(:,i) = tau_e(:)
  enddo

  ! re-defines 1D models with crustal modification if necessary
  if (CRUSTAL) then
    if (REFERENCE_1D_MODEL == REFERENCE_MODEL_AK135F_NO_MUD) then
      ! redefines 1D model with crustal modification
      call define_model_ak135(CRUSTAL)
    else if (REFERENCE_1D_MODEL == REFERENCE_MODEL_1066A) then
      ! redefines 1D model with crustal modification
      call define_model_1066a(CRUSTAL)
    else if (REFERENCE_1D_MODEL == REFERENCE_MODEL_1DREF) then
      ! redefines 1D model with crustal modification
      call define_model_1dref(CRUSTAL)
    else if (REFERENCE_1D_MODEL == REFERENCE_MODEL_SEA1D) then
      ! redefines 1D model with crustal modification
      call define_model_sea1d(CRUSTAL)
    endif
  endif

  end subroutine model_attenuation_setup

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_attenuation_getstored_tau(Qmu_in, T_c_source, tau_s, tau_e, AM_V, AM_S, AS_V)

! includes min_period, max_period, and N_SLS

  use constants

  implicit none

! model_attenuation_variables
  type model_attenuation_variables
    sequence
    double precision min_period, max_period
    double precision                          :: QT_c_source        ! Source Frequency
    double precision, dimension(:), pointer   :: Qtau_s             ! tau_sigma
    double precision, dimension(:), pointer   :: QrDisc             ! Discontinuities Defined
    double precision, dimension(:), pointer   :: Qr                 ! Radius
    double precision, dimension(:), pointer   :: Qmu                ! Shear Attenuation
    double precision, dimension(:,:), pointer :: Qtau_e             ! tau_epsilon
    double precision, dimension(:), pointer   :: Qomsb, Qomsb2      ! one_minus_sum_beta
    double precision, dimension(:,:), pointer :: Qfc, Qfc2          ! factor_common
    double precision, dimension(:), pointer   :: Qsf, Qsf2          ! scale_factor
    integer, dimension(:), pointer            :: Qrmin              ! Max and Mins of idoubling
    integer, dimension(:), pointer            :: Qrmax              ! Max and Mins of idoubling
    integer, dimension(:), pointer            :: interval_Q                 ! Steps
    integer                                   :: Qn                 ! Number of points
    integer dummy_pad ! padding 4 bytes to align the structure
  end type model_attenuation_variables

  type (model_attenuation_variables) AM_V
! model_attenuation_variables

! model_attenuation_storage_var
  type model_attenuation_storage_var
    sequence
    double precision, dimension(:,:), pointer :: tau_e_storage
    double precision, dimension(:), pointer :: Qmu_storage
    integer Q_resolution
    integer Q_max
  end type model_attenuation_storage_var

  type (model_attenuation_storage_var) AM_S
! model_attenuation_storage_var

! attenuation_simplex_variables
  type attenuation_simplex_variables
    sequence
    double precision Q  ! Q     = Desired Value of Attenuation or Q
    double precision iQ ! iQ    = 1/Q
    double precision, dimension(:), pointer ::  f
    ! f = Frequencies at which to evaluate the solution
    double precision, dimension(:), pointer :: tau_s
    ! tau_s = Tau_sigma defined by the frequency range and
    !             number of standard linear solids
    integer nf          ! nf    = Number of Frequencies
    integer nsls        ! nsls  = Number of Standard Linear Solids
  end type attenuation_simplex_variables

  type(attenuation_simplex_variables) AS_V
! attenuation_simplex_variables

  double precision :: Qmu_in, T_c_source
  double precision, dimension(N_SLS) :: tau_s, tau_e

  ! local parameters
  integer :: rw

  ! READ
  rw = 1
  call model_attenuation_storage(Qmu_in, tau_e, rw, AM_S)
  if (rw > 0) return

  call attenuation_invert_by_simplex(AM_V%min_period, AM_V%max_period, N_SLS, Qmu_in, T_c_source, tau_s, tau_e, AS_V)

  ! WRITE
  rw = -1
  call model_attenuation_storage(Qmu_in, tau_e, rw, AM_S)

  end subroutine model_attenuation_getstored_tau

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_attenuation_storage(Qmu, tau_e, rw, AM_S)

  use constants

  implicit none

! model_attenuation_storage_var
  type model_attenuation_storage_var
    sequence
    double precision, dimension(:,:), pointer :: tau_e_storage
    double precision, dimension(:), pointer :: Qmu_storage
    integer Q_resolution
    integer Q_max
  end type model_attenuation_storage_var

  type (model_attenuation_storage_var) AM_S
! model_attenuation_storage_var

  double precision :: Qmu
  double precision, dimension(N_SLS) :: tau_e
  integer :: rw

  ! local parameters
  double precision :: Qmu_new
  integer :: myrank
  integer :: Qtmp
  integer, save :: first_time_called = 1

  if (first_time_called == 1) then
     first_time_called       = 0
     AM_S%Q_resolution = 10**ATTENUATION_COMP_RESOLUTION
     AM_S%Q_max        = ATTENUATION_COMP_MAXIMUM
     Qtmp         = AM_S%Q_resolution * AM_S%Q_max
     allocate(AM_S%tau_e_storage(N_SLS, Qtmp))
     allocate(AM_S%Qmu_storage(Qtmp))
     AM_S%Qmu_storage(:) = -1
  endif

  if (Qmu < 0.0d0 .OR. Qmu > AM_S%Q_max) then
    write(IMAIN,*) 'Error attenuation_storage()'
    write(IMAIN,*) 'Attenuation Value out of Range: ', Qmu
    write(IMAIN,*) 'Attenuation Value out of Range: Min, Max ', 0, AM_S%Q_max
    call flush_IMAIN()
    ! stop
    call world_rank(myrank)
    call exit_MPI(myrank, 'Attenuation Value out of Range')
  endif

  if (rw > 0 .AND. Qmu == 0.0d0) then
     Qmu = 0.0d0;
     tau_e(:) = 0.0d0;
     return
  endif
  ! Generate index for Storage Array
  ! and Recast Qmu using this index

  ! by default: resolution is Q_resolution = 10
  ! converts Qmu to an array integer index:
  ! e.g. Qmu = 150.31 -> Qtmp = 150.31 * 10 = int( 1503.10 ) = 1503
  Qtmp    = int(Qmu * dble(AM_S%Q_resolution))

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

  subroutine attenuation_source_frequency(omega_not, min_period, max_period)
  ! Determine the Source Frequency

  implicit none

  double precision omega_not
  double precision f1, f2
  double precision min_period, max_period

  f1 = 1.0d0 / max_period
  f2 = 1.0d0 / min_period

  omega_not =  1.0e+03 * 10.0d0**(0.5 * (log10(f1) + log10(f2)))

  end subroutine attenuation_source_frequency

!
!-------------------------------------------------------------------------------------------------
!

  subroutine attenuation_tau_sigma(tau_s, n, min_period, max_period)
  ! Set the Tau_sigma (tau_s) to be equally spaced in log10 frequency

  use constants

  implicit none

  integer n
  double precision tau_s(n)
  double precision min_period, max_period
  double precision f1, f2
  double precision exp1, exp2
  double precision dexpval
  integer i

  f1 = 1.0d0 / max_period
  f2 = 1.0d0 / min_period

  exp1 = log10(f1)
  exp2 = log10(f2)

  dexpval = (exp2-exp1) / ((n*1.0d0) - 1)
  do i = 1,n
    tau_s(i) = 1.0 / (PI * 2.0d0 * 10**(exp1 + (i - 1)* 1.0d0 *dexpval))
  enddo

  end subroutine attenuation_tau_sigma

!
!-------------------------------------------------------------------------------------------------
!

  subroutine attenuation_invert_by_simplex(t2, t1, n, Q_real, omega_not, tau_s, tau_e, AS_V)

  use constants

  implicit none

! attenuation_simplex_variables
  type attenuation_simplex_variables
    sequence
    double precision Q  ! Q     = Desired Value of Attenuation or Q
    double precision iQ ! iQ    = 1/Q
    double precision, dimension(:), pointer ::  f
    ! f = Frequencies at which to evaluate the solution
    double precision, dimension(:), pointer :: tau_s
    ! tau_s = Tau_sigma defined by the frequency range and
    !             number of standard linear solids
    integer nf          ! nf    = Number of Frequencies
    integer nsls        ! nsls  = Number of Standard Linear Solids
  end type attenuation_simplex_variables

  type(attenuation_simplex_variables) AS_V
! attenuation_simplex_variables

  ! Input / Output
  integer myrank
  double precision  t1, t2
  double precision  Q_real
  double precision  omega_not
  integer  n
  double precision, dimension(n)   :: tau_s, tau_e

  ! Internal
  integer i, iterations, err,prnt
  double precision f1, f2, exp1,exp2,dexpval, min_value
  double precision, allocatable, dimension(:) :: f
  integer, parameter :: nf = 100
  double precision, external :: attenuation_eval

  ! Values to be passed into the simplex minimization routine
  iterations = -1
  min_value  = -1.0e-4
  err        = 0
  prnt       = 0

  allocate(f(nf))
  ! Determine the min and max frequencies
  f1 = 1.0d0 / t1
  f2 = 1.0d0 / t2

  ! Determine the exponents of the frequencies
  exp1 = log10(f1)
  exp2 = log10(f2)

  if (f2 < f1 .OR. Q_real < 0.0d0 .OR. n < 1) then
    call world_rank(myrank)
    call exit_MPI(myrank, 'frequencies flipped or Q less than zero or N_SLS < 0')
  endif

  ! Determine the Source frequency
  omega_not =  1.0e+03 * 10.0d0**(0.5 * (log10(f1) + log10(f2)))

  ! Determine the Frequencies at which to compare solutions
  !   The frequencies should be equally spaced in log10 frequency
  do i = 1,nf
     f(i) = exp1 + ((i-1)*1.0d0 * (exp2-exp1) / ((nf-1)*1.0d0))
  enddo

  ! Set the Tau_sigma (tau_s) to be equally spaced in log10 frequency
  dexpval = (exp2-exp1) / ((n*1.0d0) - 1)
  do i = 1,n
     tau_s(i) = 1.0 / (PI * 2.0d0 * 10**(exp1 + (i - 1)* 1.0d0 *dexpval))
  enddo

  ! Shove the parameters into the module
  call attenuation_simplex_setup(nf,n,f,Q_real,tau_s,AS_V)

  ! Set the Tau_epsilon (tau_e) to an initial value at omega*tau = 1
  ! tan_delta = 1/Q = (tau_e - tau_s)/(2 * sqrt(tau e*tau_s))
  !    if we assume tau_e =~ tau_s
  !    we get the equation below
  do i = 1,n
     tau_e(i) = tau_s(i) + (tau_s(i) * 2.0d0/Q_real)
  enddo

  ! Run a simplex search to determine the optimum values of tau_e
  call fminsearch(attenuation_eval, tau_e, n, iterations, min_value, prnt, err,AS_V)
  if (err > 0) then
     write(*,*)'Search did not converge for an attenuation of ', Q_real
     write(*,*)'    Iterations: ', iterations
     write(*,*)'    Min Value:  ', min_value
     write(*,*)'    Aborting program'
     call world_rank(myrank)
     call exit_MPI(myrank,'attenuation_simplex: Search for Strain relaxation times did not converge')
  endif
  deallocate(f)

  call attenuation_simplex_finish(AS_V)

  end subroutine attenuation_invert_by_simplex

!
!-------------------------------------------------------------------------------------------------
!

  subroutine attenuation_simplex_finish(AS_V)

  implicit none

! attenuation_simplex_variables
  type attenuation_simplex_variables
    sequence
    double precision Q  ! Q     = Desired Value of Attenuation or Q
    double precision iQ ! iQ    = 1/Q
    double precision, dimension(:), pointer ::  f
    ! f = Frequencies at which to evaluate the solution
    double precision, dimension(:), pointer :: tau_s
    ! tau_s = Tau_sigma defined by the frequency range and
    !             number of standard linear solids
    integer nf          ! nf    = Number of Frequencies
    integer nsls        ! nsls  = Number of Standard Linear Solids
  end type attenuation_simplex_variables

  type(attenuation_simplex_variables) AS_V
! attenuation_simplex_variables

  deallocate(AS_V%f)
  deallocate(AS_V%tau_s)

end subroutine attenuation_simplex_finish

!   - Inserts necessary parameters into the module attenuation_simplex_variables
!   - See module for explanation
subroutine attenuation_simplex_setup(nf_in,nsls_in,f_in,Q_in,tau_s_in,AS_V)

  implicit none

! attenuation_simplex_variables
  type attenuation_simplex_variables
    sequence
    double precision Q  ! Q     = Desired Value of Attenuation or Q
    double precision iQ ! iQ    = 1/Q
    double precision, dimension(:), pointer ::  f
    ! f = Frequencies at which to evaluate the solution
    double precision, dimension(:), pointer :: tau_s
    ! tau_s = Tau_sigma defined by the frequency range and
    !             number of standard linear solids
    integer nf          ! nf    = Number of Frequencies
    integer nsls        ! nsls  = Number of Standard Linear Solids
  end type attenuation_simplex_variables

  type(attenuation_simplex_variables) AS_V
! attenuation_simplex_variables

  integer nf_in, nsls_in
  double precision Q_in
  double precision, dimension(nf_in)   :: f_in
  double precision, dimension(nsls_in) :: tau_s_in

  allocate(AS_V%f(nf_in))
  allocate(AS_V%tau_s(nsls_in))

  AS_V%nf    = nf_in
  AS_V%nsls  = nsls_in
  AS_V%f     = f_in
  AS_V%Q     = Q_in
  AS_V%iQ    = 1.0d0/AS_V%Q
  AS_V%tau_s = tau_s_in

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
  integer nf, nsls
  double precision, dimension(nf)   :: f
  double precision, dimension(nsls) :: tau_s, tau_e
  ! Output
  double precision, dimension(nf)   :: A,B

  integer i,j
  double precision w, demon

  A(:) = 1.0d0 -  nsls*1.0d0
  B(:) = 0.0d0
  do i = 1,nf
     w = 2.0d0 * PI * 10**f(i)
     do j = 1,nsls
        demon = 1.0d0 + w**2 * tau_s(j)**2
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
  double precision function attenuation_eval(Xin,AS_V)

  implicit none

! attenuation_simplex_variables
  type attenuation_simplex_variables
    sequence
    double precision Q  ! Q     = Desired Value of Attenuation or Q
    double precision iQ ! iQ    = 1/Q
    double precision, dimension(:), pointer ::  f
    ! f = Frequencies at which to evaluate the solution
    double precision, dimension(:), pointer :: tau_s
    ! tau_s = Tau_sigma defined by the frequency range and
    !             number of standard linear solids
    integer nf          ! nf    = Number of Frequencies
    integer nsls        ! nsls  = Number of Standard Linear Solids
  end type attenuation_simplex_variables

  type(attenuation_simplex_variables) AS_V
! attenuation_simplex_variables

   ! Input
  double precision, dimension(AS_V%nsls) :: Xin
  double precision, dimension(AS_V%nsls) :: tau_e

  double precision, dimension(AS_V%nf)   :: A, B, tan_delta

  integer i
  double precision xi, iQ2

  tau_e = Xin

  call attenuation_maxwell(AS_V%nf,AS_V%nsls,AS_V%f,AS_V%tau_s,tau_e,B,A)

  tan_delta = B / A

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
!                 1 => Function Evaluation exceeded limit
!                 2 => Iterations exceeded limit
!
!     See Matlab fminsearch
  subroutine fminsearch(funk, x, n, itercount, tolf, prnt, err, AS_V)

  implicit none

! attenuation_simplex_variables
  type attenuation_simplex_variables
    sequence
    double precision Q  ! Q     = Desired Value of Attenuation or Q
    double precision iQ ! iQ    = 1/Q
    double precision, dimension(:), pointer ::  f
    ! f = Frequencies at which to evaluate the solution
    double precision, dimension(:), pointer :: tau_s
    ! tau_s = Tau_sigma defined by the frequency range and
    !             number of standard linear solids
    integer nf          ! nf    = Number of Frequencies
    integer nsls        ! nsls  = Number of Standard Linear Solids
  end type attenuation_simplex_variables

  type(attenuation_simplex_variables) AS_V
! attenuation_simplex_variables

  ! Input
  double precision, external :: funk

  integer n
  double precision x(n) ! Also Output
  integer itercount, prnt, err
  double precision tolf

  !Internal
  integer i,j, how
  integer, parameter :: none             = 0
  integer, parameter :: initial          = 1
  integer, parameter :: expand           = 2
  integer, parameter :: reflect          = 3
  integer, parameter :: contract_outside = 4
  integer, parameter :: contract_inside  = 5
  integer, parameter :: shrink           = 6

  integer maxiter, maxfun
  integer func_evals
  double precision tolx

  double precision rho, chi, psi, sigma
  double precision xin(n), y(n), v(n,n+1), fv(n+1)
  double precision vtmp(n,n+1)
  double precision usual_delta, zero_term_delta
  double precision xbar(n), xr(n), fxr, xe(n), fxe, xc(n), fxc, fxcc, xcc(n)
  integer place(n+1)

  double precision max_size_simplex, max_value

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

  fv(1) = funk(xin,AS_V)

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
     fv(j+1) = funk(x,AS_V)
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

  do while (func_evals < maxfun .AND. itercount < maxiter)

     if (max_size_simplex(v,n) <= tolx .AND. &
          max_value(fv,n+1) <= tolf) then
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
     fxr = funk(x,AS_V)
     func_evals = func_evals + 1
     if (fxr < fv(1)) then
        ! Calculate the expansion point
        xe = (1 + rho*chi)*xbar - rho*chi*v(:,n+1)
        x = xe
        fxe = funk(x,AS_V)
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
              fxc = funk(x,AS_V)
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
              fxcc = funk(x,AS_V)
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
                 fv(j) = funk(x,AS_V)
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
!      Returns:
!         Xi = max( || fv(1)- fv(i) || ) for i=2:n
!
  double precision function max_value(fv,n)
  implicit none
  integer n
  double precision fv(n)

  integer i
  double precision m, z

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
!     Returns:
!       Xi = max( max( || v(:,1) - v(:,i) || ) ) for i=2:n+1
!
  double precision function max_size_simplex(v,n)
  implicit none
  integer n
  double precision v(n,n+1)

  integer i,j
  double precision m, z

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

