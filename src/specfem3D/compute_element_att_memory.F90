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

! we switch between vectorized and non-vectorized version by using pre-processor flag FORCE_VECTORIZATION
! and macros INDEX_IJK, DO_LOOP_IJK, ENDDO_LOOP_IJK defined in config.fh
#include "config.fh"



!--------------------------------------------------------------------------------------------
!
! crust/mantle region
!
!--------------------------------------------------------------------------------------------


  subroutine compute_element_att_memory_cm(ispec,R_xx,R_yy,R_xy,R_xz,R_yz, &
                                           vx,vy,vz,vnspec,factor_common, &
                                           alphaval,betaval,gammaval, &
                                           c44store,muvstore, &
                                           epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz, &
                                           epsilondev_loc)
! crust mantle
! update memory variables based upon the Runge-Kutta scheme

! convention for attenuation
! term in xx = 1
! term in yy = 2
! term in xy = 3
! term in xz = 4
! term in yz = 5
! term in zz not computed since zero trace
! This is because we only implement Q_\mu attenuation and not Q_\kappa.
! Note that this does *NOT* imply that there is no attenuation for P waves
! because for Q_\kappa = infinity one gets (see for instance Dahlen and Tromp (1998)
! equation (9.59) page 350): Q_\alpha = Q_\mu * 3 * (V_p/V_s)^2 / 4
! therefore Q_\alpha is not zero; for instance for V_p / V_s = sqrt(3)
! we get Q_\alpha = (9 / 4) * Q_\mu = 2.25 * Q_\mu

  use constants_solver

  implicit none

  ! element id
  integer :: ispec

  ! attenuation
  ! memory variables for attenuation
  ! memory variables R_ij are stored at the local rather than global level
  ! to allow for optimization of cache access by compiler
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_ATTENUATION) :: R_xx,R_yy,R_xy,R_xz,R_yz

  ! variable sized array variables
  integer :: vx,vy,vz,vnspec

  real(kind=CUSTOM_REAL), dimension(vx,vy,vz,N_SLS,vnspec) :: factor_common
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: alphaval,betaval,gammaval

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE) :: c44store
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ISO_MANTLE) :: muvstore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: &
    epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,5) :: epsilondev_loc

! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: factor_common_c44_muv
  integer :: i_SLS

#ifdef FORCE_VECTORIZATION
  integer :: ijk
#else
  integer :: i,j,k
#endif

  ! use Runge-Kutta scheme to march in time

  ! get coefficients for that standard linear solid
  ! IMPROVE we use mu_v here even if there is some anisotropy
  ! IMPROVE we should probably use an average value instead


  do i_SLS = 1,N_SLS

    ! reformatted R_memory to handle large factor_common and reduced [alpha,beta,gamma]val
    if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then

      if (ANISOTROPIC_3D_MANTLE_VAL) then

        DO_LOOP_IJK
          factor_common_c44_muv(INDEX_IJK) = factor_common(INDEX_IJK,i_SLS,ispec) * c44store(INDEX_IJK,ispec)
        ENDDO_LOOP_IJK

      else

        DO_LOOP_IJK
          factor_common_c44_muv(INDEX_IJK) = factor_common(INDEX_IJK,i_SLS,ispec) * muvstore(INDEX_IJK,ispec)
        ENDDO_LOOP_IJK

      endif

    else

      if (ANISOTROPIC_3D_MANTLE_VAL) then

        DO_LOOP_IJK
          factor_common_c44_muv(INDEX_IJK) = factor_common(1,1,1,i_SLS,ispec) * c44store(INDEX_IJK,ispec)
        ENDDO_LOOP_IJK

      else

        DO_LOOP_IJK
          factor_common_c44_muv(INDEX_IJK) = factor_common(1,1,1,i_SLS,ispec) * muvstore(INDEX_IJK,ispec)
        ENDDO_LOOP_IJK

      endif
    endif

    ! updates memory variables
    DO_LOOP_IJK

      R_xx(INDEX_IJK,i_SLS,ispec) = alphaval(i_SLS) * R_xx(INDEX_IJK,i_SLS,ispec) + factor_common_c44_muv(INDEX_IJK) * &
          (betaval(i_SLS) * epsilondev_xx(INDEX_IJK,ispec) + gammaval(i_SLS) * epsilondev_loc(INDEX_IJK,1))

      R_yy(INDEX_IJK,i_SLS,ispec) = alphaval(i_SLS) * R_yy(INDEX_IJK,i_SLS,ispec) + factor_common_c44_muv(INDEX_IJK) * &
          (betaval(i_SLS) * epsilondev_yy(INDEX_IJK,ispec) + gammaval(i_SLS) * epsilondev_loc(INDEX_IJK,2))

      R_xy(INDEX_IJK,i_SLS,ispec) = alphaval(i_SLS) * R_xy(INDEX_IJK,i_SLS,ispec) + factor_common_c44_muv(INDEX_IJK) * &
          (betaval(i_SLS) * epsilondev_xy(INDEX_IJK,ispec) + gammaval(i_SLS) * epsilondev_loc(INDEX_IJK,3))

      R_xz(INDEX_IJK,i_SLS,ispec) = alphaval(i_SLS) * R_xz(INDEX_IJK,i_SLS,ispec) + factor_common_c44_muv(INDEX_IJK) * &
          (betaval(i_SLS) * epsilondev_xz(INDEX_IJK,ispec) + gammaval(i_SLS) * epsilondev_loc(INDEX_IJK,4))

      R_yz(INDEX_IJK,i_SLS,ispec) = alphaval(i_SLS) * R_yz(INDEX_IJK,i_SLS,ispec) + factor_common_c44_muv(INDEX_IJK) * &
          (betaval(i_SLS) * epsilondev_yz(INDEX_IJK,ispec) + gammaval(i_SLS) * epsilondev_loc(INDEX_IJK,5))

    ENDDO_LOOP_IJK

  enddo ! i_SLS

  end subroutine compute_element_att_memory_cm


!
!--------------------------------------------------------------------------------------------
!

  subroutine compute_element_att_memory_cm_lddrk(ispec,R_xx,R_yy,R_xy,R_xz,R_yz, &
                                                 R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk, &
                                                 vx,vy,vz,vnspec,factor_common, &
                                                 c44store,muvstore, &
                                                 epsilondev_loc, &
                                                 deltat)
! crust mantle
! update memory variables based upon the Runge-Kutta scheme

! convention for attenuation
! term in xx = 1
! term in yy = 2
! term in xy = 3
! term in xz = 4
! term in yz = 5
! term in zz not computed since zero trace
! This is because we only implement Q_\mu attenuation and not Q_\kappa.
! Note that this does *NOT* imply that there is no attenuation for P waves
! because for Q_\kappa = infinity one gets (see for instance Dahlen and Tromp (1998)
! equation (9.59) page 350): Q_\alpha = Q_\mu * 3 * (V_p/V_s)^2 / 4
! therefore Q_\alpha is not zero; for instance for V_p / V_s = sqrt(3)
! we get Q_\alpha = (9 / 4) * Q_\mu = 2.25 * Q_\mu

  use constants_solver
  use specfem_par,only: tau_sigma_CUSTOM_REAL,istage

  implicit none

  ! element id
  integer :: ispec

  ! attenuation
  ! memory variables for attenuation
  ! memory variables R_ij are stored at the local rather than global level
  ! to allow for optimization of cache access by compiler
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_ATTENUATION) :: R_xx,R_yy,R_xy,R_xz,R_yz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_ATTENUATION) :: &
    R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk

  ! variable sized array variables
  integer :: vx,vy,vz,vnspec

  real(kind=CUSTOM_REAL), dimension(vx,vy,vz,N_SLS,vnspec) :: factor_common

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE) :: c44store
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ISO_MANTLE) :: muvstore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,5) :: epsilondev_loc

  real(kind=CUSTOM_REAL) :: deltat

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: factor_common_c44_muv
  integer :: i_SLS

#ifdef FORCE_VECTORIZATION
  integer :: ijk
#else
  integer :: i,j,k
#endif

  ! use Runge-Kutta scheme to march in time

  ! get coefficients for that standard linear solid
  ! IMPROVE we use mu_v here even if there is some anisotropy
  ! IMPROVE we should probably use an average value instead

  do i_SLS = 1,N_SLS

    ! reformatted R_memory to handle large factor_common and reduced [alpha,beta,gamma]val
    if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
      if (ANISOTROPIC_3D_MANTLE_VAL) then

        DO_LOOP_IJK
          factor_common_c44_muv(INDEX_IJK) = factor_common(INDEX_IJK,i_SLS,ispec) * c44store(INDEX_IJK,ispec)
        ENDDO_LOOP_IJK
      else
        DO_LOOP_IJK
          factor_common_c44_muv(INDEX_IJK) = factor_common(INDEX_IJK,i_SLS,ispec) * muvstore(INDEX_IJK,ispec)
        ENDDO_LOOP_IJK
      endif

    else

      if (ANISOTROPIC_3D_MANTLE_VAL) then
        DO_LOOP_IJK
          factor_common_c44_muv(INDEX_IJK) = factor_common(1,1,1,i_SLS,ispec) * c44store(INDEX_IJK,ispec)
        ENDDO_LOOP_IJK
      else
        DO_LOOP_IJK
          factor_common_c44_muv(INDEX_IJK) = factor_common(1,1,1,i_SLS,ispec) * muvstore(INDEX_IJK,ispec)
        ENDDO_LOOP_IJK
      endif

    endif

    ! updates memory variables
    DO_LOOP_IJK

      R_xx_lddrk(INDEX_IJK,i_SLS,ispec) = ALPHA_LDDRK(istage) * R_xx_lddrk(INDEX_IJK,i_SLS,ispec) &
                       + deltat * ( factor_common_c44_muv(INDEX_IJK) * epsilondev_loc(INDEX_IJK,1) &
                                    - R_xx(INDEX_IJK,i_SLS,ispec)*(1._CUSTOM_REAL/tau_sigma_CUSTOM_REAL(i_SLS)) )

      R_yy_lddrk(INDEX_IJK,i_SLS,ispec) = ALPHA_LDDRK(istage) * R_yy_lddrk(INDEX_IJK,i_SLS,ispec) &
                       + deltat * ( factor_common_c44_muv(INDEX_IJK) * epsilondev_loc(INDEX_IJK,2) &
                                    - R_yy(INDEX_IJK,i_SLS,ispec)*(1._CUSTOM_REAL/tau_sigma_CUSTOM_REAL(i_SLS)) )

      R_xy_lddrk(INDEX_IJK,i_SLS,ispec) = ALPHA_LDDRK(istage) * R_xy_lddrk(INDEX_IJK,i_SLS,ispec) &
                       + deltat * ( factor_common_c44_muv(INDEX_IJK) * epsilondev_loc(INDEX_IJK,3) &
                                    - R_xy(INDEX_IJK,i_SLS,ispec)*(1._CUSTOM_REAL/tau_sigma_CUSTOM_REAL(i_SLS)) )

      R_xz_lddrk(INDEX_IJK,i_SLS,ispec) = ALPHA_LDDRK(istage) * R_xz_lddrk(INDEX_IJK,i_SLS,ispec) &
                       + deltat * ( factor_common_c44_muv(INDEX_IJK) * epsilondev_loc(INDEX_IJK,4) &
                                    - R_xz(INDEX_IJK,i_SLS,ispec)*(1._CUSTOM_REAL/tau_sigma_CUSTOM_REAL(i_SLS)) )

      R_yz_lddrk(INDEX_IJK,i_SLS,ispec) = ALPHA_LDDRK(istage) * R_yz_lddrk(INDEX_IJK,i_SLS,ispec) &
                       + deltat * ( factor_common_c44_muv(INDEX_IJK) * epsilondev_loc(INDEX_IJK,5) &
                                    - R_yz(INDEX_IJK,i_SLS,ispec)*(1._CUSTOM_REAL/tau_sigma_CUSTOM_REAL(i_SLS)) )

      R_xx(INDEX_IJK,i_SLS,ispec) = R_xx(INDEX_IJK,i_SLS,ispec) + BETA_LDDRK(istage) * R_xx_lddrk(INDEX_IJK,i_SLS,ispec)
      R_yy(INDEX_IJK,i_SLS,ispec) = R_yy(INDEX_IJK,i_SLS,ispec) + BETA_LDDRK(istage) * R_yy_lddrk(INDEX_IJK,i_SLS,ispec)
      R_xy(INDEX_IJK,i_SLS,ispec) = R_xy(INDEX_IJK,i_SLS,ispec) + BETA_LDDRK(istage) * R_xy_lddrk(INDEX_IJK,i_SLS,ispec)
      R_xz(INDEX_IJK,i_SLS,ispec) = R_xz(INDEX_IJK,i_SLS,ispec) + BETA_LDDRK(istage) * R_xz_lddrk(INDEX_IJK,i_SLS,ispec)
      R_yz(INDEX_IJK,i_SLS,ispec) = R_yz(INDEX_IJK,i_SLS,ispec) + BETA_LDDRK(istage) * R_yz_lddrk(INDEX_IJK,i_SLS,ispec)

    ENDDO_LOOP_IJK

  enddo ! i_SLS

  end subroutine compute_element_att_memory_cm_lddrk



!--------------------------------------------------------------------------------------------
!
! inner core region
!
!--------------------------------------------------------------------------------------------


  subroutine compute_element_att_memory_ic(ispec,R_xx,R_yy,R_xy,R_xz,R_yz, &
                                           vx,vy,vz,vnspec,factor_common, &
                                           alphaval,betaval,gammaval, &
                                           muvstore, &
                                           epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz, &
                                           epsilondev_loc)
! inner core
! update memory variables based upon the Runge-Kutta scheme

! convention for attenuation
! term in xx = 1
! term in yy = 2
! term in xy = 3
! term in xz = 4
! term in yz = 5
! term in zz not computed since zero trace
! This is because we only implement Q_\mu attenuation and not Q_\kappa.
! Note that this does *NOT* imply that there is no attenuation for P waves
! because for Q_\kappa = infinity one gets (see for instance Dahlen and Tromp (1998)
! equation (9.59) page 350): Q_\alpha = Q_\mu * 3 * (V_p/V_s)^2 / 4
! therefore Q_\alpha is not zero; for instance for V_p / V_s = sqrt(3)
! we get Q_\alpha = (9 / 4) * Q_\mu = 2.25 * Q_\mu

  use constants_solver

  implicit none

  ! element id
  integer :: ispec

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_ATTENUATION) :: R_xx,R_yy,R_xy,R_xz,R_yz

  ! variable sized array variables
  integer :: vx,vy,vz,vnspec

  real(kind=CUSTOM_REAL), dimension(vx,vy,vz,N_SLS,vnspec) :: factor_common
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: alphaval,betaval,gammaval

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: muvstore

!  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: epsilondev
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: &
    epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,5) :: epsilondev_loc

! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: factor_common_use

  integer :: i_SLS

#ifdef FORCE_VECTORIZATION
  integer :: ijk
#else
  integer :: i,j,k
#endif

  ! use Runge-Kutta scheme to march in time

  ! get coefficients for that standard linear solid
  ! IMPROVE we use mu_v here even if there is some anisotropy
  ! IMPROVE we should probably use an average value instead

  ! note: epsilondev_loc is calculated based on displ( n + 1 ), thus corresponds to strain at time (n + 1)
  !       epsilondev_xx,.. are stored from previous step, thus corresponds now to strain at time n

  do i_SLS = 1,N_SLS

    ! reformatted R_memory to handle large factor_common and reduced [alpha,beta,gamma]val
    if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then

      DO_LOOP_IJK
        factor_common_use(INDEX_IJK) = factor_common(INDEX_IJK,i_SLS,ispec) * muvstore(INDEX_IJK,ispec)
      ENDDO_LOOP_IJK

    else

      DO_LOOP_IJK
        factor_common_use(INDEX_IJK) = factor_common(1,1,1,i_SLS,ispec) * muvstore(INDEX_IJK,ispec)
      ENDDO_LOOP_IJK

    endif

    ! updates memory variables
    DO_LOOP_IJK

      R_xx(INDEX_IJK,i_SLS,ispec) = alphaval(i_SLS) * R_xx(INDEX_IJK,i_SLS,ispec) + factor_common_use(INDEX_IJK) * &
           (betaval(i_SLS) * epsilondev_xx(INDEX_IJK,ispec) + gammaval(i_SLS) * epsilondev_loc(INDEX_IJK,1))

      R_yy(INDEX_IJK,i_SLS,ispec) = alphaval(i_SLS) * R_yy(INDEX_IJK,i_SLS,ispec) + factor_common_use(INDEX_IJK) * &
           (betaval(i_SLS) * epsilondev_yy(INDEX_IJK,ispec) + gammaval(i_SLS) * epsilondev_loc(INDEX_IJK,2))

      R_xy(INDEX_IJK,i_SLS,ispec) = alphaval(i_SLS) * R_xy(INDEX_IJK,i_SLS,ispec) + factor_common_use(INDEX_IJK) * &
           (betaval(i_SLS) * epsilondev_xy(INDEX_IJK,ispec) + gammaval(i_SLS) * epsilondev_loc(INDEX_IJK,3))

      R_xz(INDEX_IJK,i_SLS,ispec) = alphaval(i_SLS) * R_xz(INDEX_IJK,i_SLS,ispec) + factor_common_use(INDEX_IJK) * &
           (betaval(i_SLS) * epsilondev_xz(INDEX_IJK,ispec) + gammaval(i_SLS) * epsilondev_loc(INDEX_IJK,4))

      R_yz(INDEX_IJK,i_SLS,ispec) = alphaval(i_SLS) * R_yz(INDEX_IJK,i_SLS,ispec) + factor_common_use(INDEX_IJK) * &
           (betaval(i_SLS) * epsilondev_yz(INDEX_IJK,ispec) + gammaval(i_SLS) * epsilondev_loc(INDEX_IJK,5))

    ENDDO_LOOP_IJK

  enddo ! N_SLS

  end subroutine compute_element_att_memory_ic


!
!--------------------------------------------------------------------------------------------
!

  subroutine compute_element_att_memory_ic_lddrk(ispec,R_xx,R_yy,R_xy,R_xz,R_yz, &
                                                 R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk, &
                                                 vx,vy,vz,vnspec,factor_common, &
                                                 muvstore, &
                                                 epsilondev_loc, &
                                                 deltat)
! inner core
! update memory variables based upon the Runge-Kutta scheme

! convention for attenuation
! term in xx = 1
! term in yy = 2
! term in xy = 3
! term in xz = 4
! term in yz = 5
! term in zz not computed since zero trace
! This is because we only implement Q_\mu attenuation and not Q_\kappa.
! Note that this does *NOT* imply that there is no attenuation for P waves
! because for Q_\kappa = infinity one gets (see for instance Dahlen and Tromp (1998)
! equation (9.59) page 350): Q_\alpha = Q_\mu * 3 * (V_p/V_s)^2 / 4
! therefore Q_\alpha is not zero; for instance for V_p / V_s = sqrt(3)
! we get Q_\alpha = (9 / 4) * Q_\mu = 2.25 * Q_\mu

  use constants_solver
  use specfem_par,only: tau_sigma_CUSTOM_REAL,istage

  implicit none

  ! element id
  integer :: ispec

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_ATTENUATION) :: &
    R_xx,R_yy,R_xy,R_xz,R_yz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_ATTENUATION) :: &
    R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk

  ! variable sized array variables
  integer :: vx,vy,vz,vnspec

  real(kind=CUSTOM_REAL), dimension(vx,vy,vz,N_SLS,vnspec) :: factor_common

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: muvstore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,5) :: epsilondev_loc

  real(kind=CUSTOM_REAL) :: deltat

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: factor_common_use

  integer :: i_SLS

#ifdef FORCE_VECTORIZATION
  integer :: ijk
#else
  integer :: i,j,k
#endif

  ! use Runge-Kutta scheme to march in time

  ! get coefficients for that standard linear solid
  ! IMPROVE we use mu_v here even if there is some anisotropy
  ! IMPROVE we should probably use an average value instead

  ! note: epsilondev_loc is calculated based on displ( n + 1 ), thus corresponds to strain at time (n + 1)
  !       epsilondev_xx,.. are stored from previous step, thus corresponds now to strain at time n

  do i_SLS = 1,N_SLS

    ! reformatted R_memory to handle large factor_common and reduced [alpha,beta,gamma]val
    if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
      DO_LOOP_IJK
        factor_common_use(INDEX_IJK) = factor_common(INDEX_IJK,i_SLS,ispec) * muvstore(INDEX_IJK,ispec)
      ENDDO_LOOP_IJK
    else
      DO_LOOP_IJK
        factor_common_use(INDEX_IJK) = factor_common(1,1,1,i_SLS,ispec) * muvstore(INDEX_IJK,ispec)
      ENDDO_LOOP_IJK
    endif

    ! updates memory variables
    DO_LOOP_IJK

      R_xx_lddrk(INDEX_IJK,i_SLS,ispec) = ALPHA_LDDRK(istage) * R_xx_lddrk(INDEX_IJK,i_SLS,ispec) &
        + deltat * ( factor_common_use(INDEX_IJK)*epsilondev_loc(INDEX_IJK,1) &
                     - R_xx(INDEX_IJK,i_SLS,ispec)*(1.0_CUSTOM_REAL/tau_sigma_CUSTOM_REAL(i_SLS)) )

      R_yy_lddrk(INDEX_IJK,i_SLS,ispec) = ALPHA_LDDRK(istage) * R_yy_lddrk(INDEX_IJK,i_SLS,ispec) &
        + deltat * ( factor_common_use(INDEX_IJK)*epsilondev_loc(INDEX_IJK,2) &
                     - R_yy(INDEX_IJK,i_SLS,ispec)*(1.0_CUSTOM_REAL/tau_sigma_CUSTOM_REAL(i_SLS)) )

      R_xy_lddrk(INDEX_IJK,i_SLS,ispec) = ALPHA_LDDRK(istage) * R_xy_lddrk(INDEX_IJK,i_SLS,ispec) &
        + deltat * ( factor_common_use(INDEX_IJK)*epsilondev_loc(INDEX_IJK,3) &
                     - R_xy(INDEX_IJK,i_SLS,ispec)*(1.0_CUSTOM_REAL/tau_sigma_CUSTOM_REAL(i_SLS)) )

      R_xz_lddrk(INDEX_IJK,i_SLS,ispec) = ALPHA_LDDRK(istage) * R_xz_lddrk(INDEX_IJK,i_SLS,ispec) &
        + deltat * ( factor_common_use(INDEX_IJK)*epsilondev_loc(INDEX_IJK,4) &
                     - R_xz(INDEX_IJK,i_SLS,ispec)*(1.0_CUSTOM_REAL/tau_sigma_CUSTOM_REAL(i_SLS)) )

      R_yz_lddrk(INDEX_IJK,i_SLS,ispec) = ALPHA_LDDRK(istage) * R_yz_lddrk(INDEX_IJK,i_SLS,ispec) &
        + deltat * ( factor_common_use(INDEX_IJK)*epsilondev_loc(INDEX_IJK,5) &
                     - R_yz(INDEX_IJK,i_SLS,ispec)*(1.0_CUSTOM_REAL/tau_sigma_CUSTOM_REAL(i_SLS)) )

      R_xx(INDEX_IJK,i_SLS,ispec) = R_xx(INDEX_IJK,i_SLS,ispec) + BETA_LDDRK(istage) * R_xx_lddrk(INDEX_IJK,i_SLS,ispec)
      R_yy(INDEX_IJK,i_SLS,ispec) = R_yy(INDEX_IJK,i_SLS,ispec) + BETA_LDDRK(istage) * R_yy_lddrk(INDEX_IJK,i_SLS,ispec)
      R_xy(INDEX_IJK,i_SLS,ispec) = R_xy(INDEX_IJK,i_SLS,ispec) + BETA_LDDRK(istage) * R_xy_lddrk(INDEX_IJK,i_SLS,ispec)
      R_xz(INDEX_IJK,i_SLS,ispec) = R_xz(INDEX_IJK,i_SLS,ispec) + BETA_LDDRK(istage) * R_xz_lddrk(INDEX_IJK,i_SLS,ispec)
      R_yz(INDEX_IJK,i_SLS,ispec) = R_yz(INDEX_IJK,i_SLS,ispec) + BETA_LDDRK(istage) * R_yz_lddrk(INDEX_IJK,i_SLS,ispec)

    ENDDO_LOOP_IJK

  enddo

  end subroutine compute_element_att_memory_ic_lddrk


!
!--------------------------------------------------------------------------------------------
!
! helper functions
!
!
!daniel debug: att - debug update
!
!  subroutine compute_element_att_mem_up_cm(ispec,i,j,k, &
!                                              R_xx_loc,R_yy_loc,R_xy_loc,R_xz_loc,R_yz_loc, &
!                                              epsilondev_loc,c44_muv,is_backward_field)
!! crust mantle
!! update memory variables based upon the Runge-Kutta scheme
!
!
!!daniel: att - debug update
!  use specfem_par,only: tau_sigma_dble,deltat,b_deltat
!
!  use specfem_par_crustmantle,only: factor_common=>factor_common_crust_mantle
!
!  use constants_solver
!
!  implicit none
!
!  ! element id
!  integer :: ispec,i,j,k
!
!  ! attenuation
!  ! memory variables for attenuation
!  ! memory variables R_ij are stored at the local rather than global level
!  ! to allow for optimization of cache access by compiler
!!  real(kind=CUSTOM_REAL), dimension(5,N_SLS) :: R_memory_loc
!  real(kind=CUSTOM_REAL), dimension(N_SLS) :: R_xx_loc
!  real(kind=CUSTOM_REAL), dimension(N_SLS) :: R_yy_loc
!  real(kind=CUSTOM_REAL), dimension(N_SLS) :: R_xy_loc
!  real(kind=CUSTOM_REAL), dimension(N_SLS) :: R_xz_loc
!  real(kind=CUSTOM_REAL), dimension(N_SLS) :: R_yz_loc
!
!  real(kind=CUSTOM_REAL), dimension(5) :: epsilondev_loc
!  real(kind=CUSTOM_REAL) :: c44_muv
!
!  logical :: is_backward_field
!  double precision :: dt,kappa
!
!! local parameters
!  real(kind=CUSTOM_REAL) :: factor_common_c44_muv
!  integer :: i_SLS
!
!  if (.not. is_backward_field) then
!    dt = dble(deltat)
!  else
!    ! backward/reconstruction: reverse time
!    dt = dble(b_deltat)
!  endif
!
!  do i_SLS = 1,N_SLS
!
!    ! Runge-Kutta scheme to update memory variables R(t)
!    if (.false.) then
!! classical RK 4:       R'(t) =  - 1/tau * R(t)
!!
!! Butcher RK4:
!! 0     |
!! 1/2   | 1/2
!! 1/2   | 0    1/2
!! 1     | 0          1
!! -----------------------------------------------------------------------------
!!         1/6  1/3   1/3   1/6
!    kappa = - dt/tau_sigma_dble(i_SLS)
!
!    R_xx_loc(i_SLS) = R_xx_loc(i_SLS) * &
!      (1.0d0 + kappa*(1.d0 + 0.5d0*kappa*(1.d0 + 1.0d0/6.0d0*kappa*(1.d0 + 1.0d0/24.0d0*kappa))))
!    R_yy_loc(i_SLS) = R_yy_loc(i_SLS) * &
!      (1.0d0 + kappa*(1.d0 + 0.5d0*kappa*(1.d0 + 1.0d0/6.0d0*kappa*(1.d0 + 1.0d0/24.0d0*kappa))))
!    R_xy_loc(i_SLS) = R_xy_loc(i_SLS) * &
!      (1.0d0 + kappa*(1.d0 + 0.5d0*kappa*(1.d0 + 1.0d0/6.0d0*kappa*(1.d0 + 1.0d0/24.0d0*kappa))))
!    R_xz_loc(i_SLS) = R_xz_loc(i_SLS) * &
!      (1.0d0 + kappa*(1.d0 + 0.5d0*kappa*(1.d0 + 1.0d0/6.0d0*kappa*(1.d0 + 1.0d0/24.0d0*kappa))))
!    R_yz_loc(i_SLS) = R_yz_loc(i_SLS) * &
!      (1.0d0 + kappa*(1.d0 + 0.5d0*kappa*(1.d0 + 1.0d0/6.0d0*kappa*(1.d0 + 1.0d0/24.0d0*kappa))))
!    endif
!
!    ! reformatted R_memory to handle large factor_common and reduced [alpha,beta,gamma]val
!    if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
!      factor_common_c44_muv = factor_common(i,j,k,i_SLS,ispec) * c44_muv
!    else
!      factor_common_c44_muv = factor_common(1,1,1,i_SLS,ispec) * c44_muv
!    endif
!
!    ! adds contributions from current strain
!    R_xx_loc(i_SLS) = R_xx_loc(i_SLS) + 0.5d0 * dt * dble(factor_common_c44_muv) * dble(epsilondev_loc(1))
!    R_yy_loc(i_SLS) = R_yy_loc(i_SLS) + 0.5d0 * dt * dble(factor_common_c44_muv) * dble(epsilondev_loc(2))
!    R_xy_loc(i_SLS) = R_xy_loc(i_SLS) + 0.5d0 * dt * dble(factor_common_c44_muv) * dble(epsilondev_loc(3))
!    R_xz_loc(i_SLS) = R_xz_loc(i_SLS) + 0.5d0 * dt * dble(factor_common_c44_muv) * dble(epsilondev_loc(4))
!    R_yz_loc(i_SLS) = R_yz_loc(i_SLS) + 0.5d0 * dt * dble(factor_common_c44_muv) * dble(epsilondev_loc(5))
!
!  enddo ! i_SLS
!
!  end subroutine compute_element_att_mem_up_cm
!

