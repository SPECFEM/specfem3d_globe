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

!-------------------------------
!
! 1D Berkeley model
!
! 1D reference model for SEMUCB - WM1
! radially anisotropic shear-wave model
!
!-------------------------------

module model_1dberkeley_par

  ! Added by < FM> Feb. 2022
  use constants, only: A3d_folder

  implicit none

  ! number of layers in model1Dberkeley.dat
  integer :: NR_REF_BERKELEY
  integer :: NR_inner_core_berk
  integer :: NR_outer_core_berk
  integer :: NR_water_berk
  integer :: ifanis_berk
  integer :: tref_berk
  integer :: ifdeck_berk

  ! model_1dberkeley_variables
  double precision, dimension(:), allocatable :: &
    Mref_V_radius_berkeley, &
    Mref_V_density_berkeley, &
    Mref_V_vpv_berkeley, &
    Mref_V_vph_berkeley, &
    Mref_V_vsv_berkeley, &
    Mref_V_vsh_berkeley, &
    Mref_V_eta_berkeley, &
    Mref_V_Qkappa_berkeley, &
    Mref_V_Qmu_berkeley

  ! Utpal Kumar, Feb, 2022
  ! define the berkeley 1D model
  character (len=100) :: berk_model1D = trim(A3d_folder) // 'model1D.dat'

  integer :: modemohoberk = -1

end module model_1dberkeley_par

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_1dberkeley_broadcast()

! reads and broadcasts berkeley 1D model

  use model_1dberkeley_par
  use constants, only: myrank

  implicit none

  ! define the berkeley 1D model
  ! Utpal Kumar, Feb, 2022

  ! main process reads in model
  if (myrank == 0) call read_1dberkeley()

  ! broadcast header values
  call BCAST_ALL_SINGLEI(ifanis_berk       )
  call BCAST_ALL_SINGLEI(tref_berk         )
  call BCAST_ALL_SINGLEI(ifdeck_berk       )
  call BCAST_ALL_SINGLEI(NR_REF_BERKELEY   )
  call BCAST_ALL_SINGLEI(NR_inner_core_berk)
  call BCAST_ALL_SINGLEI(NR_outer_core_berk)
  call BCAST_ALL_SINGLEI(NR_water_berk     )

  ! allocate arrays
  if (.not. allocated(Mref_V_radius_berkeley)) then
    allocate(Mref_V_radius_berkeley(NR_REF_BERKELEY), &
             Mref_V_density_berkeley(NR_REF_BERKELEY), &
             Mref_V_vpv_berkeley(NR_REF_BERKELEY), &
             Mref_V_vsv_berkeley(NR_REF_BERKELEY), &
             Mref_V_Qkappa_berkeley(NR_REF_BERKELEY), &
             Mref_V_Qmu_berkeley(NR_REF_BERKELEY), &
             Mref_V_vph_berkeley(NR_REF_BERKELEY), &
             Mref_V_vsh_berkeley(NR_REF_BERKELEY), &
             Mref_V_eta_berkeley(NR_REF_BERKELEY)      )
  endif

  ! broadcast data
  call BCAST_ALL_DP(Mref_V_radius_berkeley ,NR_REF_BERKELEY)
  call BCAST_ALL_DP(Mref_V_density_berkeley,NR_REF_BERKELEY)
  call BCAST_ALL_DP(Mref_V_vpv_berkeley    ,NR_REF_BERKELEY)
  call BCAST_ALL_DP(Mref_V_vph_berkeley    ,NR_REF_BERKELEY)
  call BCAST_ALL_DP(Mref_V_vsv_berkeley    ,NR_REF_BERKELEY)
  call BCAST_ALL_DP(Mref_V_vsh_berkeley    ,NR_REF_BERKELEY)
  call BCAST_ALL_DP(Mref_V_eta_berkeley    ,NR_REF_BERKELEY)
  call BCAST_ALL_DP(Mref_V_Qkappa_berkeley ,NR_REF_BERKELEY)
  call BCAST_ALL_DP(Mref_V_Qmu_berkeley    ,NR_REF_BERKELEY)

  !! Utpal Kumar, Feb, 2022
  ! determine node index where crust starts
  if (myrank == 0) then
    call determine_1dberkeley_moho_node(modemohoberk)
    !debug
    !print *,'debug: [model_1dberkeley_broadcast] modemohoberk = ',modemohoberk
  endif
  ! broadcasts to all other processes
  call bcast_all_singlei(modemohoberk)

  end subroutine model_1dberkeley_broadcast

!
!--------------------------------------------------------------------------------------------------
!

  subroutine read_1dberkeley()

  use model_1dberkeley_par

  implicit none

  integer :: i,ier
  integer, parameter :: lunit = 54
  character (len=100) :: title

  ! only main process reads in file

  open(lunit,file=trim(berk_model1D),status='old',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening file: ',trim(berk_model1D)
    stop 'Error opening file model1D.dat for berkeley reference model'
  endif

  read(lunit,100,iostat=ier) title
  read(lunit,*  ,iostat=ier) ifanis_berk, tref_berk, ifdeck_berk

  read(lunit,*  ,iostat=ier) NR_REF_BERKELEY, NR_inner_core_berk, &
                             NR_outer_core_berk, NR_water_berk

  ! allocate arrays (only on main process)
  allocate(Mref_V_radius_berkeley(NR_REF_BERKELEY), &
           Mref_V_density_berkeley(NR_REF_BERKELEY), &
           Mref_V_vpv_berkeley(NR_REF_BERKELEY), &
           Mref_V_vsv_berkeley(NR_REF_BERKELEY), &
           Mref_V_Qkappa_berkeley(NR_REF_BERKELEY), &
           Mref_V_Qmu_berkeley(NR_REF_BERKELEY), &
           Mref_V_vph_berkeley(NR_REF_BERKELEY), &
           Mref_V_vsh_berkeley(NR_REF_BERKELEY), &
           Mref_V_eta_berkeley(NR_REF_BERKELEY)      )

  ! reads data
  do i = 1,NR_REF_BERKELEY
    read(lunit,*) Mref_V_radius_berkeley(i), &
                  Mref_V_density_berkeley(i), &
                  Mref_V_vpv_berkeley(i), &
                  Mref_V_vsv_berkeley(i), &
                  Mref_V_Qkappa_berkeley(i), &
                  Mref_V_Qmu_berkeley(i), &
                  Mref_V_vph_berkeley(i), &
                  Mref_V_vsh_berkeley(i), &
                  Mref_V_eta_berkeley(i)
  enddo

  close(lunit)

  !
  ! reading formats
  !
100 format(a80)
!105 format(f8.0, 3f9.2, 2f9.1, 2f9.2, f9.5)

  end subroutine read_1dberkeley

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_1dberkeley(x,rho,vpv,vph,vsv,vsh,eta,Qkappa,Qmu,iregion_code,CRUSTAL)

  use model_1dberkeley_par
  use constants, only: PI,GRAV,EARTH_RHOAV,EARTH_R,EARTH_R_KM, &
    IREGION_INNER_CORE,IREGION_OUTER_CORE,IREGION_CRUST_MANTLE

  implicit none

! model_1dref_variables

! input:
! dimensionless radius x

! output: non-dimensionalized
!
! mass density             : rho
! compressional wave speed : vpv
! compressional wave speed : vph
! shear wave speed         : vsv
! shear wave speed         : vsh
! dimensionless parameter  : eta
! shear quality factor     : Qmu
! bulk quality factor      : Qkappa

  double precision,intent(in) :: x
  double precision,intent(inout) :: rho,vpv,vph,vsv,vsh,eta,Qmu,Qkappa
  integer,intent(in) :: iregion_code
  logical,intent(in) :: CRUSTAL

  ! local parameters
  double precision :: r,frac,scaleval
  integer :: i
  logical, parameter :: mimic_native_specfem = .true.

  ! compute real physical radius in meters
  r = x * EARTH_R

  i = 1
  do while(r >= Mref_V_radius_berkeley(i) .and. i /= NR_REF_BERKELEY)
    i = i + 1
  enddo

  ! make sure we stay in the right region
  if (mimic_native_specfem .and. iregion_code == IREGION_INNER_CORE .and. i > NR_inner_core_berk) i = NR_inner_core_berk

  if (mimic_native_specfem .and. iregion_code == IREGION_OUTER_CORE .and. i < NR_inner_core_berk+2) i = NR_inner_core_berk+2
  if (mimic_native_specfem .and. iregion_code == IREGION_OUTER_CORE .and. i > NR_outer_core_berk) i = NR_outer_core_berk

  if (mimic_native_specfem .and. iregion_code == IREGION_CRUST_MANTLE .and. i < NR_outer_core_berk+2) i = NR_outer_core_berk+2

  ! if crustal model is used, mantle gets expanded up to surface
  ! for any depth less than 24.4 km, values from mantle below moho are taken
  if (mimic_native_specfem .and. CRUSTAL .and. i > modemohoberk) then
    i = modemohoberk ! Warning : may need to be changed if file is modified !
  endif

  if (i == 1) then
    ! first layer in inner core
    rho    = Mref_V_density_berkeley(i)
    vpv    = Mref_V_vpv_berkeley(i)
    vph    = Mref_V_vph_berkeley(i)
    vsv    = Mref_V_vsv_berkeley(i)
    vsh    = Mref_V_vsh_berkeley(i)
    eta    = Mref_V_eta_berkeley(i)
    Qkappa = Mref_V_Qkappa_berkeley(i)
    Qmu    = Mref_V_Qmu_berkeley(i)
  else
    ! interpolates between one layer below to actual radius layer,
    ! that is from radius_ref(i-1) to r using the values at i-1 and i
    frac = (r-Mref_V_radius_berkeley(i-1))/(Mref_V_radius_berkeley(i)-Mref_V_radius_berkeley(i-1))
    ! interpolated model parameters
    rho = Mref_V_density_berkeley(i-1)  + frac * (Mref_V_density_berkeley(i)- Mref_V_density_berkeley(i-1))
    vpv = Mref_V_vpv_berkeley(i-1)      + frac * (Mref_V_vpv_berkeley(i)    - Mref_V_vpv_berkeley(i-1)    )
    vph = Mref_V_vph_berkeley(i-1)      + frac * (Mref_V_vph_berkeley(i)    - Mref_V_vph_berkeley(i-1)    )
    vsv = Mref_V_vsv_berkeley(i-1)      + frac * (Mref_V_vsv_berkeley(i)    - Mref_V_vsv_berkeley(i-1)    )
    vsh = Mref_V_vsh_berkeley(i-1)      + frac * (Mref_V_vsh_berkeley(i)    - Mref_V_vsh_berkeley(i-1)    )
    eta = Mref_V_eta_berkeley(i-1)      + frac * (Mref_V_eta_berkeley(i)    - Mref_V_eta_berkeley(i-1)    )
    Qkappa = Mref_V_Qkappa_berkeley(i-1)+ frac * (Mref_V_Qkappa_berkeley(i) - Mref_V_Qkappa_berkeley(i-1) )
    Qmu = Mref_V_Qmu_berkeley(i-1)      + frac * (Mref_V_Qmu_berkeley(i)    - Mref_V_Qmu_berkeley(i-1)    )
  endif

  ! make sure Vs is zero in the outer core even if roundoff errors on depth
  ! also set fictitious attenuation to a very high value (attenuation is not used in the fluid)
  if (mimic_native_specfem .and. iregion_code == IREGION_OUTER_CORE) then
    vsv = 0.d0
    vsh = 0.d0
    Qkappa = 3000.d0
    Qmu = 3000.d0
  endif

  ! non-dimensionalize
  ! time scaling (s^{-1}) is done with scaleval
  scaleval = dsqrt(PI*GRAV*EARTH_RHOAV)

  rho = rho / EARTH_RHOAV
  vpv = vpv / (EARTH_R * scaleval)
  vph = vph / (EARTH_R * scaleval)
  vsv = vsv / (EARTH_R * scaleval)
  vsh = vsh / (EARTH_R * scaleval)

  end subroutine model_1dberkeley

!
!--------------------------------------------------------------------------------------------------
!

!! Utpal Kumar, Feb, 2022
!! subroutine to decide whether the moho node has already been computed

  subroutine est_moho_node(estmohonode)

  use model_1dberkeley_par, only: modemohoberk

  implicit none
  integer :: estmohonode

  if (modemohoberk < 0) then
    call determine_1dberkeley_moho_node(estmohonode)
    ! print *,"Determining Moho node ",estmohonode
  else
    estmohonode = modemohoberk
  endif

  end subroutine est_moho_node

!
!--------------------------------------------------------------------------------------------------
!

!! Utpal Kumar, Feb 2022
!! subroutine to get the moho node number

  subroutine determine_1dberkeley_moho_node(mohonodebrk)

  use model_1dberkeley_par, only: berk_model1D

  implicit none

  integer, intent(out) :: mohonodebrk

  ! local parameter
  integer, parameter :: FID = 10
  character (len=100) :: CTMP

  double precision, allocatable :: radius(:),density(:),vpv(:),vsv(:),qkappa(:),qmu(:),vph(:),vsh(:),eta(:)
  double precision, allocatable :: derivdensity(:), derivVp(:), derivVs(:)
  double precision :: maxderivdensity

  double precision, parameter :: tol = 2.0d0**(-5)
  integer, parameter :: totalheadersval = 3

  integer :: i,IERR
  integer :: num_lines, totalDiscont
  integer :: j

  ! initializes moho node
  mohonodebrk = 1

  ! model_file = berk_model1D
  num_lines = 0 !initialize num of nodes in the file
  open(FID, file=trim(berk_model1D), status="old",iostat=IERR)
  if (IERR /= 0) then
    print *,'Error opening file: ',trim(berk_model1D)
    stop 'Error opening berkeley 1d model'
  endif

  ! Get number of lines
  do i = 1, totalheadersval
    read(FID,*)   !! skip the header
  enddo

  do while (IERR == 0)
    num_lines = num_lines + 1
    read(FID,*,iostat=IERR) CTMP
  enddo
  num_lines = num_lines - 1

  ! Allocate array of strings
  allocate(radius(num_lines), density(num_lines), vpv(num_lines), vsv(num_lines), qkappa(num_lines))
  allocate(qmu(num_lines), vph(num_lines),vsh(num_lines),eta(num_lines))
  allocate(derivdensity(num_lines), derivVp(num_lines), derivVs(num_lines))

  ! Read the file content
  rewind(FID)
  do i = 1, totalheadersval
    read(FID,*)   !! skip the header
  enddo
  do i = 1, num_lines
    read(FID,*) radius(i), density(i), vpv(i), vsv(i), qkappa(i),qmu(i), vph(i),vsh(i),eta(i)  !Read the data
  enddo
  close(FID)

  ! find the discontinuities
  totalDiscont = 0
  do i = 1, num_lines-1
    if (abs(radius(i+1)-radius(i)) < tol) then
      derivdensity(i) = density(i+1)-density(i)
      derivVp(i) = vpv(i+1)-vpv(i)
      derivVs(i) = vsv(i+1)-vsv(i)

      totalDiscont = totalDiscont + 1
    else
      derivdensity(i) = (density(i+1)-density(i))/(radius(i+1)-radius(i))
      derivVp(i) = (vpv(i+1)-vpv(i))/(radius(i+1)-radius(i))
      derivVs(i) = (vsv(i+1)-vsv(i))/(radius(i+1)-radius(i))
    endif
  enddo

  ! Determine the Mohorovicic discontinuity node
  ! Conditions to select the moho discontinuity:
  ! 1. Radius don't change, hence discontinuity
  ! 2. Vsv(i) and Vsv(i+1) > 0
  ! 3. delta Vsv > 0
  ! 4. max depth of 90km
  ! 5. max density change within 90 km from surface
  maxderivdensity = 0.d0

  j = 1
  do i = 1, num_lines-1
    if (abs(radius(i+1)-radius(i)) < tol) then
      if ((abs(vsv(i)) > tol) .and. (abs(vsv(i+1)) > tol) .and. (abs(derivVs(i)) > tol) &
          .and. (abs(radius(i)-6371000.) < 90000.)) then

        if (abs(derivdensity(i)) > maxderivdensity) then
          maxderivdensity = abs(derivdensity(i))
          mohonodebrk = i
        endif
      endif
      j = j + 1
    endif
  enddo

  deallocate(radius,density,vpv,vsv,qkappa,qmu,vph,vsh,eta)
  deallocate(derivdensity, derivVp, derivVs)

  end subroutine determine_1dberkeley_moho_node

!
!--------------------------------------------------------------------------------------------------
!

  !! subroutine to determine moho radius from reference 1D earth model
  subroutine determine_1dberkeley_moho_radius(moho_radius)

  use model_1dberkeley_par, only: berk_model1D

  use constants, only: EARTH_R_KM

  implicit none

  double precision, intent(out) :: moho_radius

  ! local parameters
  integer, parameter :: FID = 10
  character (len=100) :: CTMP

  double precision, allocatable :: radius(:),density(:),vpv(:),vsv(:),qkappa(:),qmu(:),vph(:),vsh(:),eta(:)
  double precision, allocatable :: derivdensity(:), derivVp(:), derivVs(:)
  double precision :: maxderivdensity
  double precision :: moho_r

  double precision, parameter :: tol = 2.0d0**(-5)
  integer, parameter :: totalheadersval = 3

  integer :: i, IERR
  integer :: num_lines, totalDiscont

  ! initializes
  moho_r = 0.d0

  ! model_file = berk_model1D

  num_lines = 0 !initialize num of nodes in the file
  open(FID, file=trim(berk_model1D), status="old",iostat=IERR)
  if (IERR /= 0) then
    print *,'Error opening file: ',trim(berk_model1D)
    stop 'Error opening berkeley 1d model'
  endif

  ! Get number of lines
  do i = 1, totalheadersval
    read(FID,*)   !! skip the header
  enddo

  do while (IERR == 0)
    num_lines = num_lines + 1
    read(FID,*,iostat=IERR) CTMP
  enddo
  num_lines = num_lines - 1

  ! Allocate array of strings
  allocate(radius(num_lines), density(num_lines), vpv(num_lines), vsv(num_lines), qkappa(num_lines))
  allocate(qmu(num_lines), vph(num_lines),vsh(num_lines),eta(num_lines))
  allocate(derivdensity(num_lines), derivVp(num_lines), derivVs(num_lines))

  ! Read the file content
  rewind(FID)
  do i = 1, totalheadersval
    read(FID,*)   !! skip the header
  enddo
  do i = 1, num_lines
    read(FID,*) radius(i), density(i), vpv(i), vsv(i), qkappa(i),qmu(i), vph(i),vsh(i),eta(i)  !Read the data
  enddo
  close(FID)

  ! find the discontinuities
  totalDiscont = 0
  do i = 1, num_lines-1
    if (abs(radius(i+1)-radius(i)) < tol) then
      derivdensity(i) = density(i+1)-density(i)
      derivVp(i) = vpv(i+1)-vpv(i)
      derivVs(i) = vsv(i+1)-vsv(i)

      totalDiscont = totalDiscont + 1
    else
      derivdensity(i) = (density(i+1)-density(i))/(radius(i+1)-radius(i))
      derivVp(i) = (vpv(i+1)-vpv(i))/(radius(i+1)-radius(i))
      derivVs(i) = (vsv(i+1)-vsv(i))/(radius(i+1)-radius(i))
    endif
  enddo

  ! Determine the Mohorovicic discontinuity node
  ! Conditions to select the moho discontinuity:
  ! 1. Radius don't change, hence discontinuity
  ! 2. Vsv(i) and Vsv(i+1) > 0
  ! 3. delta Vsv > 0
  ! 4. max depth of 90km
  ! 5. max density change within 90 km from surface
  maxderivdensity = 0.d0

  do i = 1, num_lines-1
    if (abs(radius(i+1)-radius(i)) < tol) then
      if ((abs(vsv(i)) > tol) .and. (abs(vsv(i+1)) > tol) .and. (abs(derivVs(i)) > tol) &
          .and. (abs(radius(i)-6371000.) < 90000.)) then

        if (abs(derivdensity(i)) > maxderivdensity) then
          maxderivdensity = abs(derivdensity(i))
          moho_r = radius(i+1)   ! new moho radius
        endif
      endif
    endif
  enddo

  deallocate(radius,density,vpv,vsv,qkappa,qmu,vph,vsh,eta)
  deallocate(derivdensity, derivVp, derivVs)

  ! return moho radius (in km)
  moho_radius = moho_r / 1000.d0 ! convert to km

  !debug
  !print *,'debug: [determine_1dberkeley_moho_radius] moho_radius = ',moho_radius

  end subroutine determine_1dberkeley_moho_radius

!
!--------------------------------------------------------------------------------------------------
!

  subroutine determine_1dberkeley_moho_depth(moho_depth)

  implicit none
  double precision, intent(out) :: moho_depth

  ! local parameters
  double precision :: moho_radius
  double precision :: earthradius = 6371.d0

  ! get moho radius (in km)
  call determine_1dberkeley_moho_radius(moho_radius)

  ! moho depth (in km)
  moho_depth = earthradius - moho_radius

  !debug
  !print *,'debug: [determine_1dberkeley_moho_depth] moho_radius = ',moho_radius

  end subroutine determine_1dberkeley_moho_depth
