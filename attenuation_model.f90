!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 6
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!       (c) California Institute of Technology September 2006
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

!  This portion of the SPECFEM3D Code was written by:
!  Brian Savage while at
!     California Institute of Technology
!     Department of Terrestrial Magnetism / Carnegie Institute of Washington
!     Univeristy of Rhode Island
!
!  <savage@uri.edu>.
!  <savage13@gps.caltech.edu>
!  <savage13@dtm.ciw.edu>
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
!   The methodology can be found in
!   Savage and Tromp, 2006, unpublished
!

module attenuation_model_constants
  implicit none

end module attenuation_model_constants

module attenuation_model_variables
  use attenuation_model_constants
  implicit none
  double precision min_period, max_period
  double precision                              :: QT_c_source        ! Source Frequency
  double precision, dimension(:), allocatable   :: Qtau_s             ! tau_sigma

  integer                                       :: Qn                 ! Number of points
  double precision, dimension(:), allocatable   :: Qr, Qs             ! Radius and Steps
  double precision, dimension(:), allocatable   :: Qmu                ! Shear Attenuation
  double precision, dimension(:,:), allocatable :: Qtau_e             ! tau_epsilon
  double precision, dimension(:), allocatable   :: Qomsb, Qomsb2      ! one_minus_sum_beta
  double precision, dimension(:,:), allocatable :: Qfc, Qfc2          ! factor_common
  double precision, dimension(:), allocatable   :: Qsf, Qsf2          ! scale_factor

  ! a 2 at the end stands for a second derivative
  ! Qkappa                - Compressional Attenutation, User Defined
  ! Qmu                   - Shear Attenuation, User Defined
  !
  ! tau_sigma             -                       use attenuation_tau_sigma()
  !         Dependent on min and max period
  ! tau_epislon           -                       use attenuation_invert_by_simplex()
  !         Dependent on tau_sigma and Qmu
  ! [alpha,beta,gamma]val - For memory variables  use attenuation_memory_values()
  !         Dependent on tau_simga and delta_T
  ! one_minus_sum_beta    -                       use attenuation_property_values()
  !         Dependent on tau_sigma and Qmu
  ! factor_common         -                       use attenuation_property_values()
  !         Dependent on tau_sigma and Qmu
  ! scale_factor          - Scaling the velocity  use attenuation_scale_factor()
  !         Dependent on tau_sigma and Qmu

end module attenuation_model_variables

module attenuation_model_storage
  implicit none
  integer Q_resolution, Q_max
  double precision, dimension(:,:), allocatable :: tau_e_storage
  double precision, dimension(:),   allocatable :: Qmu_storage

end module attenuation_model_storage

subroutine attenuation_lookup_value(i, r)
  implicit none
  include 'constants.h'
  integer i
  double precision r
  r = dble(i) / TABLE_ATTENUATION
end subroutine attenuation_lookup_value


 subroutine attenuation_lookup_index(i, r, theta, idoubling, ell_d80, ellipticity_val)
   implicit none
   include 'mpif.h'
   include 'constants.h'
   integer i, idoubling, ier, myrank
   logical ellipticity_val
   real(kind=CUSTOM_REAL) r, theta
   real(kind=CUSTOM_REAL) cost, p20, ell_d80
   if(r < 100.d0 / R_EARTH) r = 100.d0 / R_EARTH


   i = max(1, nint(r * TABLE_ATTENUATION))
   if(ELLIPTICITY_VAL .AND. idoubling <= IFLAG_220_MOHO) then
      ! particular case of d80 which is not honored by the mesh
      ! map ellipticity back for d80 detection
      ! ystore contains theta
      cost = cos(theta)
      p20 = 0.5 * (3.0 * cost * cost - 1.0)
      r = r * (1.0 + (2.0/3.0) * ell_d80 * p20)
   endif
   if(i > NRAD_ATTENUATION) then
      write(*,*)'index: ',i,r
      call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ier)
      call exit_MPI(myrank, 'bad index in attenuation lookup table')
   endif
 end subroutine attenuation_lookup_index


! This Subroutine is Hackish.  It could probably all be moved to an input attenuation file.
! Actually all the velocities, densities and attenuations could be moved to seperate input
! files rather than be defined within the CODE
!
! All this subroutine does is define the Attenuation vs Radius and then Compute the Attenuation
! Variables (tau_sigma and tau_epslion ( or tau_mu) )
subroutine attenuation_model_setup(REFERENCE_1D_MODEL, RICB, RCMB, R670, R220, R80)
  use attenuation_model_variables
  use model_ak135_variables
  use model_1066a_variables
  implicit none
  include 'mpif.h'
  include 'constants.h'

  integer myrank
  integer REFERENCE_1D_MODEL
  double precision RICB, RCMB, R670, R220, R80
  double precision tau_e(N_SLS)

  integer i,ier
  double precision Qb
  double precision R120

  Qb = 57287.0d0
  R120 = 6251.d3
  
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ier)
  if(myrank > 0) return

  call define_model_ak135(.FALSE.)
  call define_model_1066a(.FALSE.)

  if(REFERENCE_1D_MODEL == REFERENCE_MODEL_PREM) then
     Qn = 12
  else if(REFERENCE_1D_MODEL == REFERENCE_MODEL_IASP91) then
     Qn = 12
  else if(REFERENCE_1D_MODEL == REFERENCE_MODEL_AK135) then
     Qn = NR_AK135
  else if(REFERENCE_1D_MODEL == REFERENCE_MODEL_1066a) then
     Qn = NR_1066A
  else
     call exit_MPI(myrank, 'Reference 1D Model Not recognized')
  endif

  allocate(Qr(Qn))
  allocate(Qmu(Qn))
  allocate(Qs(Qn))
  allocate(Qtau_e(N_SLS,Qn))

  if(REFERENCE_1D_MODEL == REFERENCE_MODEL_PREM) then
     Qr(:)     = (/    0.0d0,     RICB,  RICB,  RCMB,    RCMB,    R670,    R670,   R220,    R220,    R80,     R80, R_EARTH /)
     Qmu(:)    = (/   84.6d0,   84.6d0, 0.0d0, 0.0d0, 312.0d0, 312.0d0, 143.0d0, 143.0d0, 80.0d0, 80.0d0, 600.0d0, 600.0d0 /)
  else if(REFERENCE_1D_MODEL == REFERENCE_MODEL_IASP91) then
     Qr(:)     = (/    0.0d0,     RICB,  RICB,  RCMB,    RCMB,    R670,    R670,    R220,   R220,   R120,    R120, R_EARTH /)
     Qmu(:)    = (/   84.6d0,   84.6d0, 0.0d0, 0.0d0, 312.0d0, 312.0d0, 143.0d0, 143.0d0, 80.0d0, 80.0d0, 600.0d0, 600.0d0 /)
  else if(REFERENCE_1D_MODEL == REFERENCE_MODEL_AK135) then
     Qr(:)     = radius_ak135(:)
     Qmu(:)    = Qmu_ak135(:)
  else if(REFERENCE_1D_MODEL == REFERENCE_MODEL_1066a) then
     Qr(:)     = radius_1066a(:)
     Qmu(:)    = Qmu_1066a(:)
  end if

  do i = 1, Qn
     call attenuation_conversion(Qmu(i), QT_c_source, Qtau_s, tau_e)
     Qtau_e(:,i) = tau_e(:)
  end do

end subroutine attenuation_model_setup

subroutine attenuation_save_arrays(prname, iregion_code)
  use attenuation_model_variables
  implicit none
  include 'mpif.h'
  include 'constants.h'

  integer iregion_code
  character(len=150) prname
  integer ier
  integer myrank
  integer, SAVE :: virgin = 1

  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ier)
  if(myrank == 0 .AND. iregion_code == IREGION_CRUST_MANTLE .AND. virgin == 1) then
     virgin = 0
     open(unit=27,file=prname(1:len_trim(prname))//'1D_Q.bin', status='unknown', form='unformatted')
     write(27) QT_c_source
     write(27) Qtau_s
     write(27) Qn
     write(27) Qr
     write(27) Qmu
     write(27) Qtau_e
     close(27)
  else
     return
  endif

end subroutine attenuation_save_arrays

subroutine attenuation_storage(Qmu, tau_e, rw)
  use attenuation_model_storage
  implicit none
  include 'mpif.h'
  include 'constants.h'

  integer myrank, ier
  double precision Qmu
  double precision, dimension(N_SLS) :: tau_e
  integer rw

  integer Qtmp
  integer, SAVE :: virgin = 1

  if(virgin == 1) then
     virgin       = 0
     Q_resolution = 10**ATTENUATION_COMP_RESOLUTION
     Q_max        = ATTENUATION_COMP_MAXIMUM
     Qtmp         = Q_resolution * Q_max
     allocate(tau_e_storage(N_SLS, Qtmp))
     allocate(Qmu_storage(Qtmp))
     Qmu_storage(:) = -1
  endif


  if(Qmu < 0.0d0 .OR. Qmu >= Q_max) then
     write(IMAIN,*) 'Error'
     write(IMAIN,*) 'attenuation_conversion/storage()'
     write(IMAIN,*) 'Attenuation Value out of Range: ', Qmu
     write(IMAIN,*) 'Attenuation Value out of Range: Min, Max ', 0, Q_max
     call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ier)
     call exit_MPI(myrank, 'Attenuation Value out of Range')
  endif

  if(rw > 0 .AND. Qmu == 0.0d0) then
     Qmu = 0.0d0;
     tau_e(:) = 0.0d0;
     return
  endif
  ! Generate index for Storage Array
  ! and Recast Qmu using this index
  Qtmp = Qmu * Q_resolution
  Qmu = Qtmp / Q_resolution;

  if(rw > 0) then
     ! READ
     if(Qmu_storage(Qtmp) > 0) then
        ! READ SUCCESSFUL
        tau_e(:)   = tau_e_storage(:, Qtmp)
        Qmu        = Qmu_storage(Qtmp)
        rw = 1
     else
        ! READ NOT SUCCESSFUL
        rw = -1
     endif
  else
     ! WRITE SUCCESSFUL
     tau_e_storage(:,Qtmp)    = tau_e(:)
     Qmu_storage(Qtmp)        = Qmu
     rw = 1
  endif

  return

end subroutine attenuation_storage

subroutine attenuation_conversion(Qmu_in, T_c_source, tau_s, tau_e)
  use attenuation_model_variables ! includes min_period, max_period, and N_SLS
!  use attenuation_model_storage
  implicit none
  include 'constants.h'

  double precision Qmu_in, T_c_source
  double precision, dimension(N_SLS) :: tau_s, tau_e

  integer rw

  rw = 1
  call attenuation_storage(Qmu_in, tau_e, rw)
  if(rw > 0) return

  call attenuation_invert_by_simplex(min_period, max_period, N_SLS, Qmu_in, T_c_source, tau_s, tau_e)

  rw = -1
  call attenuation_storage(Qmu_in, tau_e, rw)

end subroutine attenuation_conversion

subroutine read_attenuation_model(min, max)

  use attenuation_model_variables
  implicit none
  include 'constants.h'

  integer min, max

  min_period = min * 1.0d0
  max_period = max * 1.0d0

  allocate(Qtau_s(N_SLS))
  call attenuation_tau_sigma(Qtau_s, N_SLS, min_period, max_period)
  call attenuation_source_frequency(QT_c_source, min_period, max_period)

end subroutine read_attenuation_model

subroutine attenuation_memory_values(tau_s, deltat, alphaval,betaval,gammaval)
  use attenuation_model_variables
  implicit none
  include 'constants.h'

  double precision, dimension(N_SLS) :: tau_s, alphaval, betaval,gammaval
  real(kind=CUSTOM_REAL) deltat

  double precision, dimension(N_SLS) :: tauinv

  tauinv(:) = - 1.0 / tau_s(:)

  alphaval(:)  = 1 + deltat*tauinv(:) + deltat**2*tauinv(:)**2 / 2. + &
       deltat**3*tauinv(:)**3 / 6. + deltat**4*tauinv(:)**4 / 24.
  betaval(:)   = deltat / 2. + deltat**2*tauinv(:) / 3. + deltat**3*tauinv(:)**2 / 8. + deltat**4*tauinv(:)**3 / 24.
  gammaval(:)  = deltat / 2. + deltat**2*tauinv(:) / 6. + deltat**3*tauinv(:)**2 / 24.0

end subroutine attenuation_memory_values

subroutine attenuation_scale_factor(myrank, T_c_source, tau_mu, tau_sigma, Q_mu, scale_factor)
  use attenuation_model_variables
  implicit none
  include 'constants.h'

  integer myrank
  double precision scale_factor, Q_mu, T_c_source
  double precision, dimension(N_SLS) :: tau_mu, tau_sigma

  double precision scale_t
  double precision f_c_source, w_c_source, f_0_prem
  double precision factor_scale_mu0, factor_scale_mu
  double precision a_val, b_val
  double precision big_omega
  integer i

  scale_t = ONE/dsqrt(PI*GRAV*RHOAV)

!--- compute central angular frequency of source (non dimensionalized)
  f_c_source = ONE / T_c_source
  w_c_source = TWO_PI * f_c_source

!--- non dimensionalize PREM reference of 1 second
  f_0_prem = ONE / ( ONE / scale_t)

!--- quantity by which to scale mu_0 to get mu
  factor_scale_mu0 = ONE + TWO * log(f_c_source / f_0_prem) / (PI * Q_mu)

!--- compute a, b and Omega parameters, also compute one minus sum of betas
  a_val = ONE
  b_val = ZERO

  do i = 1,N_SLS
    a_val = a_val - w_c_source * w_c_source * tau_mu(i) * &
      (tau_mu(i) - tau_sigma(i)) / (1.d0 + w_c_source * w_c_source * tau_mu(i) * tau_mu(i))
    b_val = b_val + w_c_source * (tau_mu(i) - tau_sigma(i)) / &
      (1.d0 + w_c_source * w_c_source * tau_mu(i) * tau_mu(i))
  enddo

  big_omega = a_val*(sqrt(1.d0 + b_val*b_val/(a_val*a_val))-1.d0)

!--- quantity by which to scale mu to get mu_relaxed
  factor_scale_mu = b_val * b_val / (TWO * big_omega)

!--- total factor by which to scale mu0
  scale_factor = factor_scale_mu * factor_scale_mu0

!--- check that the correction factor is close to one
  if(scale_factor < 0.9 .or. scale_factor > 1.1) then
     write(*,*)'scale factor: ', scale_factor
     call exit_MPI(myrank,'incorrect correction factor in attenuation model')
  endif

end subroutine attenuation_scale_factor

!----

subroutine attenuation_property_values(tau_s, tau_e, factor_common, one_minus_sum_beta)

  use attenuation_model_variables
  implicit none
  include 'constants.h'

  double precision, dimension(N_SLS) :: tau_s, tau_e, beta, factor_common
  double precision  one_minus_sum_beta

  double precision, dimension(N_SLS) :: tauinv
  integer i

  tauinv(:) = -1.0d0 / tau_s(:)

  beta(:) = 1.0d0 - tau_e(:) / tau_s(:)
  one_minus_sum_beta = 1.0d0

  do i = 1,N_SLS
     one_minus_sum_beta = one_minus_sum_beta - beta(i)
  enddo

  factor_common(:) = 2.0d0 * beta(:) * tauinv(:)

end subroutine attenuation_property_values

subroutine get_attenuation_model_1D(myrank, prname, iregion_code, tau_s, one_minus_sum_beta, &
                                    factor_common, scale_factor, vn,vx,vy,vz)
  use attenuation_model_variables
  implicit none
  include 'mpif.h'
  include 'constants.h'
  integer myrank, iregion_code
  character(len=150) prname
  integer vn, vx,vy,vz
  double precision, dimension(N_SLS)              :: tau_s
  double precision, dimension(vx,vy,vz,vn)        :: scale_factor, one_minus_sum_beta
  double precision, dimension(N_SLS, vx,vy,vz,vn) :: factor_common

  integer i,j,ier,rmax
  double precision scale_t
  double precision Qp1, Qpn, radius, fctmp
  double precision, dimension(:), allocatable :: Qfctmp, Qfc2tmp

  integer, SAVE :: virgin = 1

  if(myrank == 0 .AND. iregion_code == IREGION_CRUST_MANTLE .AND. virgin == 1) then
     open(unit=27, file=prname(1:len_trim(prname))//'1D_Q.bin', status='unknown', form='unformatted')
     read(27) QT_c_source
     read(27) tau_s
     read(27) Qn
     allocate(Qr(Qn))
     allocate(Qmu(Qn))
     allocate(Qtau_e(N_SLS,Qn))
     read(27) Qr
     read(27) Qmu
     read(27) Qtau_e
     close(27)
  endif

  ! Synch up after the Read
  call MPI_BCAST(QT_c_source, 1,     MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
  call MPI_BCAST(tau_s,       N_SLS, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
  call MPI_BCAST(Qn,          1,     MPI_INTEGER,          0, MPI_COMM_WORLD, ier)
  call MPI_BCAST(QT_c_source, 1,     MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
  if(myrank .NE. 0) then
     allocate(Qr(Qn))
     allocate(Qmu(Qn))
     allocate(Qtau_e(N_SLS,Qn))
  endif
  call MPI_BCAST(Qr,     Qn,       MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
  call MPI_BCAST(Qmu,    Qn,       MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
  call MPI_BCAST(Qtau_e, Qn*N_SLS, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)

  scale_t = ONE/dsqrt(PI*GRAV*RHOAV)

  ! Scale the Attenuation Values
  tau_s(:)    = tau_s(:)    / scale_t
  Qtau_e(:,:) = Qtau_e(:,:) / scale_t
  QT_c_source = 1000.0d0 / QT_c_source / scale_t
  Qr(:)       = Qr(:) / R_EARTH

  allocate(Qsf(Qn))
  allocate(Qomsb(Qn))
  allocate(Qfc(N_SLS,Qn))

  allocate(Qsf2(Qn))
  allocate(Qomsb2(Qn))
  allocate(Qfc2(N_SLS,Qn))

  allocate(Qs(Qn))

  allocate(Qfctmp(Qn))
  allocate(Qfc2tmp(Qn))

  do i = 1,Qn
     if(Qmu(i) == 0.0d0) then
        Qomsb(i) = 0.0d0
        Qfc(:,i) = 0.0d0
        Qsf(i)   = 0.0d0
     else
        call attenuation_property_values(tau_s, Qtau_e(:,i), Qfc(:,i), Qomsb(i))
        call attenuation_scale_factor(myrank, QT_c_source, Qtau_e(:,i), tau_s, Qmu(i), Qsf(i))
     endif
  enddo

  ! Determine the Spline Coefficients or Second Derivatives
  call pspline(Qr, Qsf,   Qn, Qp1, Qpn, Qsf2,   Qs)
  call pspline(Qr, Qomsb, Qn, Qp1, Qpn, Qomsb2, Qs)
  do i = 1,N_SLS
     call pspline(Qr, Qfc(i,:), Qn, Qp1, Qpn, Qfc2(i,:),Qs)
  enddo

  radius = 0.0d0
  rmax = nint(TABLE_ATTENUATION)
  do i = 1,rmax
     call attenuation_lookup_value(i, radius)
     call psplint(Qr, Qsf,   Qsf2,   Qn, radius, scale_factor(1,1,1,i),       Qs)
     call psplint(Qr, Qomsb, Qomsb2, Qn, radius, one_minus_sum_beta(1,1,1,i), Qs)
     do j = 1,N_SLS
        Qfctmp  = Qfc(j,:)
        Qfc2tmp = Qfc2(j,:)
        call psplint(Qr, Qfctmp, Qfc2tmp, Qn, radius, fctmp, Qs)
        factor_common(j,1,1,1,i) = fctmp
     enddo
  enddo
  do i = rmax+1,NRAD_ATTENUATION
     scale_factor(1,1,1,i)       = scale_factor(1,1,1,rmax)
     one_minus_sum_beta(1,1,1,i) = one_minus_sum_beta(1,1,1,rmax)
     factor_common(1,1,1,1,i)    = factor_common(1,1,1,1,rmax)
     factor_common(2,1,1,1,i)    = factor_common(2,1,1,1,rmax)
     factor_common(3,1,1,1,i)    = factor_common(3,1,1,1,rmax)
  end do

  deallocate(Qfc2)
  deallocate(Qsf2)
  deallocate(Qomsb2)
  deallocate(Qfc)
  deallocate(Qsf)
  deallocate(Qomsb)
  deallocate(Qr)
  deallocate(Qs)
  deallocate(Qtau_e)
  deallocate(Qmu)
  deallocate(Qfctmp)
  deallocate(Qfc2tmp)

  call MPI_BARRIER(MPI_COMM_WORLD, ier)

end subroutine get_attenuation_model_1D


subroutine get_attenuation_model_3D(myrank, prname, one_minus_sum_beta, factor_common, scale_factor, tau_s, vnspec)
  use attenuation_model_variables
  implicit none
  include 'constants.h'

  integer myrank, vnspec
  character(len=150) prname
  double precision, dimension(NGLLX,NGLLY,NGLLZ,vnspec)       :: one_minus_sum_beta, scale_factor
  double precision, dimension(N_SLS,NGLLX,NGLLY,NGLLZ,vnspec) :: factor_common
  double precision, dimension(N_SLS)                          :: tau_s

  integer i,j,k,ispec

  double precision, dimension(N_SLS) :: tau_e, fc
  double precision  omsb, Q_mu, sf, T_c_source, scale_t

  ! All of the following reads use the output parameters as their temporary arrays
  ! use the filename to determine the actual contents of the read

  open(unit=27, file=prname(1:len_trim(prname))//'tau_s.bin',status='old',action='read',form='unformatted')
  read(27) tau_s
  close(27)

  open(unit=27, file=prname(1:len_trim(prname))//'T_c_source.bin',status='old',action='read',form='unformatted')
  read(27) T_c_source
  close(27)

  open(unit=27, file=prname(1:len_trim(prname))//'Q.bin',status='old',action='read',form='unformatted')
  read(27) scale_factor
  close(27)

  open(unit=27, file=prname(1:len_trim(prname))//'tau_e.bin',status='old',action='read',form='unformatted')
  read(27) factor_common
  close(27)


  scale_t = ONE/dsqrt(PI*GRAV*RHOAV)

  factor_common(:,:,:,:,:) = factor_common(:,:,:,:,:) / scale_t ! This is really tau_e, not factor_common
  tau_s(:)                 = tau_s(:) / scale_t
  T_c_source               = 1000.0d0 / T_c_source
  T_c_source               = T_c_source / scale_t

  do ispec = 1, vnspec
     do k = 1, NGLLZ
        do j = 1, NGLLY
           do i = 1, NGLLX
              tau_e(:) = factor_common(:,i,j,k,ispec)
              Q_mu     = scale_factor(i,j,k,ispec)

              ! Determine the factor_common and one_minus_sum_beta from tau_s and tau_e
              call attenuation_property_values(tau_s, tau_e, fc, omsb)

              factor_common(:,i,j,k,ispec)    = fc(:)
              one_minus_sum_beta(i,j,k,ispec) = omsb

              ! Determine the "scale_factor" from tau_s, tau_e, central source frequency, and Q
              call attenuation_scale_factor(myrank, T_c_source, tau_e, tau_s, Q_mu, sf)
              scale_factor(i,j,k,ispec) = sf
           enddo
        enddo
     enddo
  enddo
end subroutine get_attenuation_model_3D

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

subroutine attenuation_tau_sigma(tau_s, n, min_period, max_period)
  ! Set the Tau_sigma (tau_s) to be equally spaced in log10 frequency
  implicit none
  integer n
  double precision tau_s(n)
  double precision min_period, max_period
  double precision f1, f2
  double precision exp1, exp2
  double precision dexp
  integer i
  double precision, parameter :: PI = 3.14159265358979d0

  f1 = 1.0d0 / max_period
  f2 = 1.0d0 / min_period

  exp1 = log10(f1)
  exp2 = log10(f2)

  dexp = (exp2-exp1) / ((n*1.0d0) - 1)
  do i = 1,n
     tau_s(i) = 1.0 / (PI * 2.0d0 * 10**(exp1 + (i - 1)* 1.0d0 *dexp))
  enddo

end subroutine attenuation_tau_sigma

module attenuation_simplex_variables
  implicit none

  ! nf    = Number of Frequencies
  ! nsls  = Number of Standard Linear Solids
  integer  nf, nsls

  ! f     = Frequencies at which to evaluate the solution
  double precision, allocatable, dimension(:) ::  f

  ! Q     = Desired Value of Attenuation or Q
  ! iQ    = 1/Q
  double precision  Q, iQ

  ! tau_s = Tau_sigma defined by the frequency range and
  !             number of standard linear solids
  double precision, allocatable, dimension(:) :: tau_s

end module attenuation_simplex_variables

subroutine attenuation_invert_by_simplex(t2, t1, n, Q_real, omega_not, tau_s, tau_e)
  implicit none
  include 'mpif.h'
  ! Input / Output
  integer myrank, ier
  double precision  t1, t2
  double precision  Q_real
  double precision  omega_not
  integer  n
  double precision, dimension(n)   :: tau_s, tau_e

  ! Internal
  integer i, iterations, err,prnt
  double precision f1, f2, exp1,exp2,dexp, min_value
  double precision, allocatable, dimension(:) :: f
  double precision, parameter :: PI = 3.14159265358979d0
  integer, parameter :: nf = 100
  double precision attenuation_eval
  EXTERNAL attenuation_eval

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

  if(f2 < f1 .OR. Q_real < 0.0d0 .OR. n < 1) then
     call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ier)     
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
  dexp = (exp2-exp1) / ((n*1.0d0) - 1)
  do i = 1,n
     tau_s(i) = 1.0 / (PI * 2.0d0 * 10**(exp1 + (i - 1)* 1.0d0 *dexp))
  enddo

  ! Shove the paramters into the module
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
  if(err > 0) then
     write(*,*)'Search did not converge for an attenuation of ', Q_real
     write(*,*)'    Iterations: ', iterations
     write(*,*)'    Min Value:  ', min_value
     write(*,*)'    Aborting program'
     call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ier)
     call exit_MPI(myrank,'attenuation_simplex: Search for Strain relaxation times did not converge')
  endif

  deallocate(f)
  call attenuation_simplex_finish()

end subroutine attenuation_invert_by_simplex


subroutine attenuation_simplex_finish()
  use attenuation_simplex_variables
  implicit none

  deallocate(f)
  deallocate(tau_s)

end subroutine attenuation_simplex_finish

!!!!!!!
! subroutine simplex_setup
!   - Inserts necessary parameters into the module attenuation_simplex_variables
!   - See module for explaination
subroutine attenuation_simplex_setup(nf_in,nsls_in,f_in,Q_in,tau_s_in)
  use attenuation_simplex_variables
  implicit none

  integer nf_in, nsls_in
  double precision Q_in
  double precision, dimension(nf_in)   :: f_in
  double precision, dimension(nsls_in) :: tau_s_in

  allocate(f(nf_in))
  allocate(tau_s(nsls_in))

  nf    = nf_in
  nsls  = nsls_in
  f     = f_in
  Q     = Q_in
  iQ    = 1.0d0/Q
  tau_s = tau_s_in

end subroutine attenuation_simplex_setup

!!!!!!!!
! subroutine attenuation_maxwell
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
!     tau_e = Tau_epislon Strain relaxation time (see References)
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
  implicit none

  ! Input
  integer nf, nsls
  double precision, dimension(nf)   :: f
  double precision, dimension(nsls) :: tau_s, tau_e
  ! Output
  double precision, dimension(nf)   :: A,B

  integer i,j
  double precision w, pi, demon

  PI = 3.14159265358979d0

  A(:) = 1.0d0 -  nsls*1.0d0
  B(:) = 0.0d0
  do i = 1,nf
     w = 2.0d0 * PI * 10**f(i)
     do j = 1,nsls
!        write(*,*)j,tau_s(j),tau_e(j)
        demon = 1.0d0 + w**2 * tau_s(j)**2
        A(i) = A(i) + ((1.0d0 + (w**2 * tau_e(j) * tau_s(j)))/ demon)
        B(i) = B(i) + ((w * (tau_e(j) - tau_s(j))) / demon)
     end do
!     write(*,*)A(i),B(i),10**f(i)
  enddo

end subroutine attenuation_maxwell

!!!!!!!!
! subroutine attenuation_eval
!    - Computes the misfit from a set of relaxation paramters
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
!           Qt_i is the desired attenuaiton at that frequency
!
!    Uses attenuation_simplex_variables to store constant values
!
!    See atteunation_simplex_setup
!
double precision function attenuation_eval(Xin)
  use attenuation_simplex_variables
  implicit none
   ! Input
  double precision, dimension(nsls) :: Xin
  double precision, dimension(nsls) :: tau_e

  double precision, dimension(nf)   :: A, B, tan_delta

  integer i
  double precision xi, iQ2

  tau_e = Xin

  call attenuation_maxwell(nf,nsls,f,tau_s,tau_e,B,A)

  tan_delta = B / A

  attenuation_eval = 0.0d0
  iQ2 = iQ**2
  do i = 1,nf
     xi = sqrt(( ( (tan_delta(i) - iQ) ** 2 ) / iQ2 ))
     attenuation_eval = attenuation_eval + xi
  enddo

end function attenuation_eval


!!!!!!!!!!!!!1
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
!            Output: Mimimized Value
!     n    = number of variables
!     itercount = Input/Output
!                 Input:  maximum number of iterations
!                         if < 0 default is used (200 * n)
!                 Output: total number of iterations on output
!     tolf      = Input/Output
!                 Input:  minimium tolerance of the function funk(x)
!                 Output: minimium value of funk(x)(i.e. "a" solution)
!     prnt      = Input
!                 3 => report every iteration
!                 4 => report every iteration, total simplex
!     err       = Output
!                 0 => Normal exeecution, converged within desired range
!                 1 => Function Evaluation exceeded limit
!                 2 => Iterations exceeded limit
!
!     See Matlab fminsearch
!
subroutine fminsearch(funk, x, n, itercount, tolf, prnt, err)
  implicit none

  ! Input
  double precision funk
  EXTERNAL funk

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


  if(itercount > 0) then
     maxiter = itercount
  else
     maxiter = 200 * n
  endif
  itercount = 0
  maxfun  = 200 * n

  if(tolf > 0.0d0) then
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
     if(y(j) .NE. 0.0d0) then
        y(j) = (1.0d0 + usual_delta) * y(j)
     else
        y(j) = zero_term_delta
     endif
     v(:,j+1) = y
     x(:) = y
     fv(j+1) = funk(x)
  enddo

  call qsort(fv,n+1,place)

  do i = 1,n+1
     vtmp(:,i) = v(:,place(i))
  enddo
  v = vtmp

  how = initial
  itercount = 1
  func_evals = n+1
  if(prnt .EQ. 3) then
     write(*,*)'Iterations   Funk Evals   Value How'
     write(*,*)itercount, func_evals, fv(1), how
  endif
  if(prnt .EQ. 4) then
     write(*,*)'How: ',how
     write(*,*)'V: ', v
     write(*,*)'fv: ',fv
     write(*,*)'evals: ',func_evals
  endif

  do while (func_evals < maxfun .AND. itercount < maxiter)

!     if(max(max(abs(v(:,2:n+1) - v(:,1)))) .LE. tolx .AND. &
!          max(abs(fv(1) - fv(2:n+1))) .LE. tolf) then

     if(max_size_simplex(v,n) .LE. tolx .AND. &
          max_value(fv,n+1) .LE. tolf) then
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
           if (how .EQ. shrink) then
              do j=2,n+1
                 v(:,j)=v(:,1)+sigma*(v(:,j) - v(:,1))
                 x(:) = v(:,j)
                 fv(j) = funk(x)
              enddo
              func_evals = func_evals + n
           endif
        endif
     endif

     call qsort(fv,n+1,place)
     do i = 1,n+1
        vtmp(:,i) = v(:,place(i))
     enddo
     v = vtmp

     itercount = itercount + 1
     if (prnt == 3) then
        write(*,*)itercount, func_evals, fv(1), how
     elseif (prnt == 4) then
        write(*,*)
        write(*,*)'How: ',how
        write(*,*)'v: ',v
        write(*,*)'fv: ',fv
        write(*,*)'evals: ',func_evals
     endif
  enddo

  if(func_evals > maxfun) then
     write(*,*)'function evaluations exceeded prescribed limit', maxfun
     err = 1
  endif
  if(itercount > maxiter) then
     write(*,*)'iterations exceeded prescribed limit', maxiter
     err = 2
  endif

666 continue
  x = v(:,1)
  tolf = fv(1)

end subroutine fminsearch


!!!!!!!
! double precision function max_value
!    - Finds the maximim value of the difference of between the first
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
     if(z > m) then
        m = z
     endif
  enddo

  max_value = m

end function max_value

!!!!!!!!
! function max_size_simplex
!   - Determines the maximum distance between two point in a simplex
!   Input
!     v  = Input
!            Simplex Verticies
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
        if(z > m) then
           m = z
        endif
     enddo
  enddo

  max_size_simplex = m

end function max_size_simplex


!!!!!!!
! subroutine qsort
!    - Implementation of a Bubble Sort Routine
!    Input
!      X = Input/Output
!         Vector to be sorted
!         dimension(n)
!      n = Input
!         Length of X
!      I = Output
!         Sorted Indicies of vecotr X
!
!      Example:
!         X = [ 4 3 1 2 ] on Input
!         I = [ 1 2 3 4 ] Computed Internally (in order)
!
!         X = [ 1 2 3 4 ] on Output
!         I = [ 3 4 2 1 ] on Output
!
subroutine qsort(X,n,I)
  implicit none
  integer n
  double precision X(n)
  integer I(n)

  integer j,k
  double precision rtmp
  integer itmp

  do j = 1,n
     I(j) = j
  enddo

  do j = 1,n
     do k = 1,n-j
        if(X(k+1) < X(k)) then
           rtmp   = X(k)
           X(k)   = X(k+1)
           X(k+1) = rtmp

           itmp   = I(k)
           I(k)   = I(k+1)
           I(k+1) = itmp
        endif
     enddo
  enddo

end subroutine qsort

! Piecewise Continuous Splines
!   - Uses the Numerical Recipes for Spline interpolation
!   - Added Steps which describes the discontinuities
!   - Steps must be repeats in the dependent variable, X
!   - Derivates at the steps are computed using the point
!     at the derivate and the closest point within that piece
!   - A point lying directly on the discontinuity will recieve the
!     value of the first or smallest piece in terms of X
!   - Beginning and Ending points of the Function become beginning
!     and ending points of the first and last splines
!   - A Step with a value of zero is undefined
!   - Works with functions with steps or no steps


subroutine psplint(xa, ya, y2a, n, x, y, steps)
  implicit none
  integer n
  double precision xa(n),ya(n),y2a(n)
  integer steps(n)
  double precision x, y

  integer i, l, n1, n2

  do i = 1,n-1
     if(steps(i+1) == 0) return
     if(x >= xa(steps(i)) .AND. x <= xa(steps(i+1))) then
        call pspline_piece(i,n1,n2,l,n,steps)
        call splint(xa(n1), ya(n1), y2a(n1), l, x, y)
        return
     endif
  end do

  return

end subroutine psplint

subroutine pspline_piece(i,n1,n2,l,n,s)
  implicit none
  integer i, n1, n2, l, n, s(n)
  n1 = s(i)+1
  if(i == 1) n1 = s(i)
  n2 = s(i+1)
  l = n2 - n1 + 1
  return
end subroutine pspline_piece

subroutine pspline(x, y, n, yp1, ypn, y2, steps)
  implicit none
  integer n
  double precision x(n),y(n),y2(n)
  double precision yp1, ypn
  integer steps(n)

  integer i,r, l, n1,n2

  steps(:) = 0

  ! Find steps in x, defining pieces
  steps(1) = 1
  r = 2
  do i = 2,n
     if(x(i) == x(i-1)) then
        steps(r) = i-1
        r = r + 1
     endif
  end do
  steps(r) = n

  ! Run spline for each piece
  do i = 1,r-1
     call pspline_piece(i,n1,n2,l,n,steps)
     ! Determine the First Derivates at Begin/End Points
     yp1 = ( y(n1+1) - y(n1) ) / ( x(n1+1) - x(n1))
     ypn = ( y(n2) - y(n2-1) ) / ( x(n2) - x(n2-1))
     call spline(x(n1),y(n1),l,yp1,ypn,y2(n1))
  enddo

end subroutine pspline



subroutine attenuation_model_1D(myrank, iregion_attenuation, Q_mu)
  implicit none
  include "constants.h"
  integer myrank
  integer iregion_attenuation
  double precision Q_mu

  ! This is the PREM Attenuation Structure
  ! check in which region we are based upon doubling flag
  select case(iregion_attenuation)
  case(IREGION_ATTENUATION_INNER_CORE) !--- inner core, target Q_mu: 84.60
     Q_mu =        84.6000000000d0
  case(IREGION_ATTENUATION_CMB_670)    !--- CMB -> d670 (no attenuation in fluid outer core), target Q_mu = 312.
     Q_mu =       312.0000000000d0
  case(IREGION_ATTENUATION_670_220)    !--- d670 -> d220, target Q_mu: 143.
     Q_mu =       143.0000000000d0
  case(IREGION_ATTENUATION_220_80)     !--- d220 -> depth of 80 km, target Q_mu:  80.
     Q_mu =        80.0000000000d0
  case(IREGION_ATTENUATION_80_SURFACE) !--- depth of 80 km -> surface, target Q_mu: 600.
     Q_mu =       600.0000000000d0
  !--- do nothing for fluid outer core (no attenuation there)
  case default
     call exit_MPI(myrank,'wrong attenuation flag in mesh')
  end select

end subroutine attenuation_model_1D

subroutine attenuation_model(x, xlat, xlon,Qmu)

! x in the radius from 0 to 1 where 0 is the center and 1 is the surface
! xlat, xlon currently not used in this routine (which uses PREM).
! The user needs to modify this routine if he wants to use
! a particular 3D attenuation model. The current version is 1D PREM.
!

!  use attenuation_model_variables

  implicit none
  include 'constants.h'

  double precision xlat, xlon, r, x, Qmu,RICB,RCMB, &
      RTOPDDOUBLEPRIME,R600,R670,R220,R771,R400,R80, ROCEAN, &
      RMOHO, RMIDDLE_CRUST
  double precision Qkappa

! dummy lines to suppress warning about variables not used
  double precision xdummy
  xdummy = xlat + xlon
! end of dummy lines to suppress warning about variables not used

  r = x * R_EARTH

  ROCEAN = 6368000.d0
  RMIDDLE_CRUST = 6356000.d0
  RMOHO = 6346600.d0
  R80  = 6291000.d0
  R220 = 6151000.d0
  R400 = 5971000.d0
  R600 = 5771000.d0
  R670 = 5701000.d0
  R771 = 5600000.d0
  RTOPDDOUBLEPRIME = 3630000.d0
  RCMB = 3480000.d0
  RICB = 1221000.d0


! PREM
!
!--- inner core
!
  if(r >= 0.d0 .and. r <= RICB) then
     Qmu=84.6d0
     Qkappa=1327.7d0
!
!--- outer core
!
  else if(r > RICB .and. r <= RCMB) then
     Qmu=0.0d0
     Qkappa=57827.0d0
     if(RCMB - r .LT. r - RICB) then
        Qmu = 312.0d0  ! CMB
     else
        Qmu = 84.6d0   ! ICB
     endif
!
!--- D" at the base of the mantle
!
  else if(r > RCMB .and. r <= RTOPDDOUBLEPRIME) then
     Qmu=312.0d0
     Qkappa=57827.0d0
!
!--- mantle: from top of D" to d670
!
  else if(r > RTOPDDOUBLEPRIME .and. r <= R771) then
     Qmu=312.0d0
     Qkappa=57827.0d0
  else if(r > R771 .and. r <= R670) then
     Qmu=312.0d0
     Qkappa=57827.0d0
!
!--- mantle: above d670
!
  else if(r > R670 .and. r <= R600) then
     Qmu=143.0d0
     Qkappa=57827.0d0
  else if(r > R600 .and. r <= R400) then
     Qmu=143.0d0
     Qkappa=57827.0d0
  else if(r > R400 .and. r <= R220) then
     Qmu=143.0d0
     Qkappa=57827.0d0
  else if(r > R220 .and. r <= R80) then
     Qmu=80.0d0
     Qkappa=57827.0d0
  else if(r > R80) then
     Qmu=600.0d0
     Qkappa=57827.0d0
  endif

!  This will be called from get_model.f90
!  call attenuation_conversion(Qmu, T_c_source, tau_s, tau_e)

end subroutine attenuation_model


! This is not used anymore

!  This subroutine works fine for large values of Q
!     However when Q becomes small ( Q < 80 ), the methodology breaks down
!  call attenuation_invert_by_SVD(myrank, min_period, max_period, N_SLS, Qmu, T_c_source, tau_s, tau_e)

!  We are now using a much robust method of attenuation parameter determination
!    Although this method may be a bit slower, it will determine the parameters
!    when the value of Q is low (Q < 80)
!   attenuation_invert_by_simplex()


! subroutine attenuation_invert_by_SVD(myrank, t2, t1, n,Q_real,omega_not,tau_s,tau_e)
!   use attenuation_model_variables
!   implicit none
!   include 'constants.h'

!   double precision  t1, t2
!   integer  myrank, n
!   double precision  Q_real
!   double precision  omega_not
!   double precision, dimension(n)   :: tau_s, tau_e

!   integer i, j, k

!   double precision, dimension(n)   :: gradient, dadp, dbdp, dqdp, x1, x2
!   double precision, dimension(n,n) :: hessian

!   double precision a, b, demon, Q_omega, Q_ratio, F_ratio, PI2, q
!   double precision f1, f2, exp1, exp2, expo, dexp, omega, df, d2qdp2

!   gradient(:)  = 0.0d0
!   hessian(:,:) = 0.0d0
!   tau_e(:)     = 0.0d0
!   tau_s(:)     = 0.0d0

!   PI2 = 6.28318530717958d0

!   f1 = 1.0d0/t1
!   f2 = 1.0d0/t2

!   if(f2 < f1) call exit_MPI(myrank, 'max frequency is less than min frequency')

!   if(Q_real < 0.0) call exit_MPI(myrank, 'Attenuation is less than zero')

!   if(n < 1) call exit_MPI(myrank, 'Number of standard linear solids is less than one')

!   omega_not = 1.0e+3 * 10.0**(0.5 * (log10(f1) + log10(f2)))

!   exp1 = log10(f1)
!   exp2 = log10(f2)

!   dexp = (exp2 - exp1) / real(n - 1.0)

!   q = 1.0 / (real(n - 1.0) * Q_real )

!   do i = 1,n
!      expo     = exp1 + (i-1) * dexp
!      omega    = PI2 * 10.0**expo
!      tau_s(i) = 1.0 / omega
!      tau_e(i) = tau_s(i) * (1.0 + q) / (1.0 - q)
!   enddo

!   x1(:) = tau_e(:) - tau_s(:)
!   x2(:) = tau_s(:)

!   exp1 = log10(f1)
!   exp2 = log10(f2)
!   dexp = (exp2 - exp1) / 100.0

!   expo = exp1 - dexp
!   do i = 1,100
!      expo = expo + dexp
!      df       = 10.0**(expo+dexp) - 10.0**(expo)
!      omega    = PI2 * 10.0**(expo)
!      a = real(1.0 - n)
!      b = 0.0
!      do j = 1,n
!         tau_e(j) = x1(j) + x2(j)
!         tau_s(j) = x2(j)
!         demon   = 1.0 + omega**2.0 * tau_s(j)**2.0
!         a       = a + (1.0 + omega**2.0 * tau_e(j) * tau_s(j)) / demon
!         b       = b + ( omega * ( tau_e(j) - tau_s(j) ) ) / demon
!         dadp(j) = omega**2.0 * tau_s(j) / demon
!         dbdp(j) = omega / demon
!      enddo

!      Q_omega = a / b
!      Q_ratio = 1.0 / Q_omega - 1.0 / Q_real
!      F_ratio = df / (f2 - f1)
!      do j = 1,n
!         dqdp(j)     = (dbdp(j) - ( b / a ) * dadp(j)) / a
!         gradient(j) = gradient(j) + 2.0 * (Q_ratio) * dqdp(j) * F_ratio
!         do k = 1,j
!            d2qdp2   = -(dadp(j) * dbdp(k) + dbdp(j) * dadp(k) - 2.0 * (b / a) * dadp(j) * dadp(k)) / (a * a)
!            hessian(j,k) = hessian(j,k) + (2.0 * dqdp(j) * dqdp(k) + 2.0 * Q_ratio * d2qdp2) * F_ratio
!            hessian(k,j) = hessian(j,k)
!         enddo
!      enddo
!   enddo

!   call invert(x1, gradient, hessian, n)
!   tau_e(:) = x1(:) + x2(:)
!   tau_s(:) = x2(:)

! end subroutine attenuation_invert_by_SVD

! !----

! subroutine invert(x,b,A,n)

!   use attenuation_model_variables
!   implicit none
!   include 'constants.h'

!   integer n

!   double precision, dimension(n)   :: x, b
!   double precision, dimension(n,n) :: A

!   integer i, j, k
!   double precision, dimension(n)   :: W, xp
!   double precision, dimension(n,n) :: V
!   double precision, dimension(n,n) :: A_inverse

!   call svdcmp_dp(A,W,V,n)

!   do i = 1,n
!      do j = 1,n
!         V(i,j) = (1.0d0 / W(i)) * A(j,i)
!      enddo
!   enddo

!   do i = 1,n
!      do j = 1,n
!         A_inverse(i,j) = 0.0d0
!         do k = 1,n
!            A_inverse(i,j) = A_inverse(i,j) + A(i,k) * V(k,j)
!         enddo
!      enddo
!   enddo

!   do i = 1,n
!      xp(i) = x(i)
!      do j = 1, n
!         xp(i) = xp(i) - A_inverse(i,j) * b(j)
!      enddo
!      x(i) = xp(i)
!   enddo

! end subroutine invert


! FUNCTION pythag_dp(a,b)
!   use attenuation_model_variables
!   IMPLICIT NONE
!   include 'constants.h'

!   double precision, INTENT(IN) :: a,b
!   double precision :: pythag_dp
!   double precision :: absa,absb
!   absa=abs(a)
!   absb=abs(b)
!   if (absa > absb) then
!      pythag_dp=absa*sqrt(1.0d0+(absb/absa)**2)
!   else
!      if (absb == 0.0d0) then
!         pythag_dp=0.0d0
!      else
!         pythag_dp=absb*sqrt(1.0+(absa/absb)**2)
!      endif
!   endif
! END FUNCTION pythag_dp

! SUBROUTINE svdcmp_dp(a,w,v,p)

!   use attenuation_model_variables

!   IMPLICIT NONE
!   include 'constants.h'


!   integer p
!   INTEGER, PARAMETER :: DP = KIND(1.0D0)
!   double precision, DIMENSION(p,p), INTENT(INOUT) :: a
!   double precision, DIMENSION(p), INTENT(OUT) :: w
!   double precision, DIMENSION(p,p), INTENT(OUT) :: v
!   INTEGER(4) :: i,its,j,k,l,m,n,nm
!   double precision :: anorm,c,f,g,h,s,scale,x,y,z
!   double precision, DIMENSION(size(a,1)) :: tempm
!   double precision, DIMENSION(size(a,2)) :: rv1,tempn
!   double precision PYTHAG_DP

!   m=size(a,1)
!   n = size(a,2)

!   g=0.0d0
!   scale=0.0d0
!   do i=1,n
!      l=i+1
!      rv1(i)=scale*g
!      g=0.0d0
!      scale=0.0d0

!      if (i <= m) then
!         scale=sum(abs(a(i:m,i)))
!         if (scale /= 0.0d0) then
!            a(i:m,i)=a(i:m,i)/scale
!            s=dot_product(a(i:m,i),a(i:m,i))
!            f=a(i,i)
!            g=-sign(sqrt(s),f)
!            h=f*g-s
!            a(i,i)=f-g
!            tempn(l:n)=matmul(a(i:m,i),a(i:m,l:n))/h

!            a(i:m,l:n)=a(i:m,l:n)+spread(a(i:m,i),dim=2,ncopies=size(tempn(l:n))) * &
!                 spread(tempn(l:n),dim=1,ncopies=size(a(i:m,i)))
!            a(i:m,i)=scale*a(i:m,i)
!         endif
!      endif
!      w(i)=scale*g
!      g=0.0d0
!      scale=0.0d0
!      if ((i <= m) .and. (i /= n)) then
!         scale=sum(abs(a(i,l:n)))
!         if (scale /= 0.0d0) then
!            a(i,l:n)=a(i,l:n)/scale
!            s=dot_product(a(i,l:n),a(i,l:n))
!            f=a(i,l)
!            g=-sign(sqrt(s),f)
!            h=f*g-s
!            a(i,l)=f-g
!            rv1(l:n)=a(i,l:n)/h
!            tempm(l:m)=matmul(a(l:m,l:n),a(i,l:n))

!            a(l:m,l:n)=a(l:m,l:n)+spread(tempm(l:m),dim=2,ncopies=size(rv1(l:n))) * &
!                 spread(rv1(l:n),dim=1,ncopies=size(tempm(l:m)))
!            a(i,l:n)=scale*a(i,l:n)
!         endif
!      endif
!   enddo
!   anorm=maxval(abs(w)+abs(rv1))

!   do i=n,1,-1
!      if (i < n) then
!         if (g /= 0.0d0) then
!            v(l:n,i)=(a(i,l:n)/a(i,l))/g
!            tempn(l:n)=matmul(a(i,l:n),v(l:n,l:n))
!            v(l:n,l:n)=v(l:n,l:n)+spread(v(l:n,i),dim=2,ncopies=size(tempn(l:n))) * &
!                 spread(tempn(l:n), dim=1, ncopies=size(v(l:n,i)))
!         endif
!         v(i,l:n)=0.0d0
!         v(l:n,i)=0.0d0
!      endif
!      v(i,i)=1.0d0
!      g=rv1(i)
!      l=i
!   enddo
!   do i=min(m,n),1,-1
!      l=i+1
!      g=w(i)
!      a(i,l:n)=0.0d0
!      if (g /= 0.0d0) then
!         g=1.0d0/g
!         tempn(l:n)=(matmul(a(l:m,i),a(l:m,l:n))/a(i,i))*g
!         a(i:m,l:n)=a(i:m,l:n)+spread(a(i:m,i),dim=2,ncopies=size(tempn(l:n))) * &
!              spread(tempn(l:n),dim=1,ncopies=size(a(i:m,i)))
!         a(i:m,i)=a(i:m,i)*g
!      else
!         a(i:m,i)=0.0d0
!      endif
!      a(i,i)=a(i,i)+1.0d0
!   enddo
!   do k=n,1,-1
!      do its=1,30
!         do l=k,1,-1
!            nm=l-1
!            if ((abs(rv1(l))+anorm) == anorm) exit
!            if ((abs(w(nm))+anorm) == anorm) then
!               c=0.0d0
!               s=1.0d0
!               do i=l,k
!                  f=s*rv1(i)
!                  rv1(i)=c*rv1(i)
!                  if ((abs(f)+anorm) == anorm) exit
!                  g=w(i)
!                  h=pythag_dp(f,g)
!                  w(i)=h
!                  h=1.0d0/h
!                  c= (g*h)
!                  s=-(f*h)
!                  tempm(1:m)=a(1:m,nm)
!                  a(1:m,nm)=a(1:m,nm)*c+a(1:m,i)*s
!                  a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
!               enddo
!               exit
!            endif
!         enddo
!         z=w(k)
!         if (l == k) then
!            if (z < 0.0d0) then
!               w(k)=-z
!               v(1:n,k)=-v(1:n,k)
!            endif
!            exit
!         endif

!         if (its == 30) stop 'svdcmp_dp: no convergence in svdcmp'

!         x=w(l)
!         nm=k-1
!         y=w(nm)
!         g=rv1(nm)
!         h=rv1(k)
!         f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0d0*h*y)
!         g=pythag_dp(f,1.0d0)
!         f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
!         c=1.0d0
!         s=1.0d0
!         do j=l,nm
!            i=j+1
!            g=rv1(i)
!            y=w(i)
!            h=s*g
!            g=c*g
!            z=pythag_dp(f,h)
!            rv1(j)=z
!            c=f/z
!            s=h/z
!            f= (x*c)+(g*s)
!            g=-(x*s)+(g*c)
!            h=y*s
!            y=y*c
!            tempn(1:n)=v(1:n,j)
!            v(1:n,j)=v(1:n,j)*c+v(1:n,i)*s
!            v(1:n,i)=-tempn(1:n)*s+v(1:n,i)*c
!            z=pythag_dp(f,h)
!            w(j)=z
!            if (z /= 0.0d0) then
!               z=1.0d0/z
!               c=f*z
!               s=h*z
!            endif
!            f= (c*g)+(s*y)
!            x=-(s*g)+(c*y)
!            tempm(1:m)=a(1:m,j)
!            a(1:m,j)=a(1:m,j)*c+a(1:m,i)*s
!            a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
!         enddo
!         rv1(l)=0.0d0
!         rv1(k)=f
!         w(k)=x
!      enddo
!   enddo

! END SUBROUTINE svdcmp_dp


