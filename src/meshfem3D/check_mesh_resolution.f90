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


  subroutine check_mesh_resolution(iregion_code,xstore,ystore,zstore, &
                                   kappavstore,kappahstore,muvstore,muhstore,rhostore)

  use constants

  use shared_parameters, only: T_min_period,ADIOS_ENABLED

  use meshfem_par, only: &
    nspec,nglob, &
    DT,myrank,ibool,idoubling, &
    xstore_glob,ystore_glob,zstore_glob,SAVE_MESH_FILES, &
    dt_max_glob,pmax_glob

  use regions_mesh_par2, only: prname

  implicit none

  integer,intent(in) :: iregion_code

  double precision,dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: xstore,ystore,zstore

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: kappavstore,kappahstore, &
    muvstore,muhstore,rhostore

  ! local parameters
  real(kind=CUSTOM_REAL) :: dt_max,dt_max_reg
  real(kind=CUSTOM_REAL) :: pmax,pmax_reg
  real(kind=CUSTOM_REAL) :: pmax_empirical
  real(kind=CUSTOM_REAL) :: elemsize_min,elemsize_max,elemsize_min_reg,elemsize_max_reg,dx
  real(kind=CUSTOM_REAL) :: cmax,cmax_reg

  real(kind=CUSTOM_REAL) :: vpmax,vsmin,vpmax_reg,vsmin_reg
  real(kind=CUSTOM_REAL) :: avg_distance,distance_min
  double precision :: deltat_suggested
  integer :: ispec,ier

  ! minimum distance between adjacent GLL points
  double precision, dimension(15), parameter :: percent_GLL = &
    (/ 100.d0, 100.d0, 50.d0, 27.639320225002102d0, 17.267316464601141d0, 11.747233803526763d0, &
       8.4888051860716516d0, 6.4129925745196719d0, 5.0121002294269914d0, 4.0233045916770571d0, 3.2999284795970416d0, &
       2.7550363888558858d0, 2.3345076678918053d0, 2.0032477366369594d0, 1.7377036748080721d0 /) / 100._CUSTOM_REAL

  ! smallest eigenvalue
  real(kind=CUSTOM_REAL) :: eig_ratio_min,eig_ratio_max
  real(kind=CUSTOM_REAL) :: eig_ratio_min_reg,eig_ratio_max_reg

  ! for file output
  logical, parameter :: USE_VTU_FORMAT = .true.           ! using VTU binary file format as output
  logical, parameter :: DEBUG_OUTPUT_FOR_ADIOS = .false.  ! outputs maximum/minimum DT/period files when ADIOS_ENABLED is .true.

  logical :: do_output
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: val_ispec_pmax,val_ispec_dt
  character(len=MAX_STRING_LEN) :: filename
  character(len=32),parameter :: region(MAX_NUM_REGIONS) = &
    (/character(len=32) :: 'crust/mantle', 'outer core', 'inner core', &
                           'transition-to-infinite', 'infinite' /)

  !debug timing
  !double precision, external :: wtime
  !double precision :: tstart,tCPU

  ! to avoid conversion warnings
  real(kind=CUSTOM_REAL), parameter :: HUGEVAL_REAL = real(HUGEVAL,kind=CUSTOM_REAL)

  ! note: the mesh and time step check is only approximative

  !debug timing
  !tstart = wtime()

  ! safety check
  if (NGLLX < 1 .or. NGLLX > 15) stop 'Invalid NGLLX value in routine check_mesh_resolution'

  ! check if we want to output files
  if (SAVE_MESH_FILES) then
    ! by default, lets the SAVE_MESH_FILES flag determine about the file output
    do_output = .true.
    ! ADIOS case
    if (ADIOS_ENABLED) then
      if (DEBUG_OUTPUT_FOR_ADIOS) then
        ! we will output vtk/vtu files nevertheless
        do_output = .true.
      else
        ! omit the file output
        do_output = .false.
      endif
    endif
  else
    ! no output requested
    do_output = .false.
  endif

  ! temporary arrays for file output
  if (do_output) then
    allocate(val_ispec_pmax(nspec),val_ispec_dt(nspec),stat=ier)
    if (ier /= 0) stop 'Error allocating val_ispec_pmax arrays'
    val_ispec_pmax(:) = 0._CUSTOM_REAL
    val_ispec_dt(:) = 0._CUSTOM_REAL
  endif

  ! initializes global min/max values only when first called
  if (iregion_code == 1) then
    dt_max_glob = HUGEVAL_REAL
    pmax_glob = - HUGEVAL_REAL
  endif

  ! statistics for this region
  elemsize_min_reg = HUGEVAL_REAL
  elemsize_max_reg = - HUGEVAL_REAL

  eig_ratio_min_reg = HUGEVAL_REAL
  eig_ratio_max_reg = - HUGEVAL_REAL

  pmax_reg = - HUGEVAL_REAL
  dt_max_reg = HUGEVAL_REAL
  cmax_reg = - HUGEVAL_REAL

  vsmin_reg = HUGEVAL_REAL
  vpmax_reg = - HUGEVAL_REAL

! openmp mesher
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ispec,vpmax,vsmin,elemsize_min,elemsize_max, &
!$OMP eig_ratio_min,eig_ratio_max, &
!$OMP avg_distance,pmax,distance_min,dx,deltat_suggested,dt_max,cmax)
!$OMP DO
  do ispec = 1,nspec

    ! suppress fictitious elements in central cube
    if (idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle

    ! determines maximum velocities at corners of this element
    call get_vpvs_minmax(vpmax,vsmin,ispec,nspec,iregion_code,kappavstore,kappahstore,muvstore,muhstore,rhostore)

    ! sets region min/max
    vsmin_reg = min(vsmin,vsmin_reg)
    vpmax_reg = max(vpmax,vpmax_reg)

    ! computes minimum size of this grid cell
    call get_elem_minmaxsize(elemsize_min,elemsize_max,ispec,nspec,xstore,ystore,zstore)

    ! sets region min/max
    elemsize_min_reg = min(elemsize_min_reg,elemsize_min)
    elemsize_max_reg = max(elemsize_max_reg,elemsize_max)

    ! gets ratio smallest/largest eigenvalue of the jacobian matrix defined for all points in this element
    call get_eigenvalues_min_ratio(eig_ratio_min,eig_ratio_max,ispec)

    eig_ratio_min_reg = min(eig_ratio_min_reg,eig_ratio_min)
    eig_ratio_max_reg = max(eig_ratio_max_reg,eig_ratio_max)

    ! debug
    !if (eig_ratio_min_reg > eig_min_ratio) then
    !  eig_ratio_min_reg = eig_min_ratio
    !  print *,'eigen: ',eig_ratio_min_reg, &
    !          'at radius ',sngl(sqrt(xstore(3,3,3,ispec)**2 + ystore(3,3,3,ispec)**2 + zstore(3,3,3,ispec)**2)*R_PLANET_KM)
    !endif

    ! largest possible minimum period such that number of points per minimum wavelength
    ! npts = ( min(vpmin,vsmin)  * pmax ) / avg_distance  is about ~ NPTS_PER_WAVELENGTH
    !
    ! note: obviously, this estimation depends on the choice of points per wavelength
    !          which is empirical at the moment.
    !          also, keep in mind that the minimum period is just an estimation and
    !          there is no such sharp cut-off period for valid synthetics.
    !          seismograms become just more and more inaccurate for periods shorter than this estimate.
    ! average distance between GLL points within this element
    avg_distance = elemsize_min / ( NGLLX - 1)  ! since NGLLX = NGLLY = NGLLZ

    ! minimum period
    pmax = avg_distance / vsmin * NPTS_PER_WAVELENGTH

    ! debug
    !if (pmax_reg < pmax) then
    !  pmax_reg = pmax
    !  print *,'minimum period = ',sngl(pmax_reg), &
    !          sngl(256.0/max(NEX_ETA,NEX_XI) * max(ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES)/90.0 * 17.0), &
    !          'vpmax = ',vpmax,'vsmin = ',vsmin, &
    !          'at radius ',sngl(sqrt(xstore(3,3,3,ispec)**2 + ystore(3,3,3,ispec)**2 + zstore(3,3,3,ispec)**2)*R_PLANET_KM)
    !endif

    ! sets region minimum period
    pmax_reg = max(pmax_reg,pmax)

    ! computes minimum and maximum distance of neighbor GLL points in this grid cell
    distance_min = elemsize_min * real(percent_GLL(NGLLX),kind=CUSTOM_REAL)

    ! distance at skewed corner across
    ! if the angle at corner less than 60 degrees, then the distance between the second GLL points (B-C) becomes
    ! smaller than the distance from corner to second GLL point (A - B or A - C):
    !
    !                   o A
    !                  / \
    !                 /   \
    !             B  o=====o C
    !
    ! note: this condition is usually not met for the case of our global meshes
    call get_min_distance_from_second_GLL_points(dx,ispec,nspec,xstore,ystore,zstore)
    if (distance_min > dx) distance_min = dx

    ! suggested timestep
    deltat_suggested = COURANT_SUGGESTED * distance_min / vpmax

    ! cut at a significant number of digits (2 digits)
    ! example: 0.0734815 -> lpow = (2 - (-1) = 3 -> 0.0730
    call get_timestep_limit_significant_digit(deltat_suggested)

    ! maximum time step size
    dt_max = real(deltat_suggested,kind=CUSTOM_REAL)

    ! debug
    !if (dt_max_reg > dt_max) then
    !  dt_max_reg = dt_max
    !  print *,'dt_max = ',dt_max_reg,sngl(DT),sngl(fac_pow),dt_cut, &
    !          'vp = ',vpmax,'distance_min = ',distance_min,eig_ratio_min, &
    !          'at radius ',sngl(sqrt(xstore(3,3,3,ispec)**2 + ystore(3,3,3,ispec)**2 + zstore(3,3,3,ispec)**2)*R_PLANET_KM)
    !endif

    ! sets region time step in this slice
    dt_max_reg = min(dt_max_reg,dt_max)

    ! Courant number
    ! based on minimum GLL point distance and maximum velocity
    ! i.e. on the maximum ratio of ( velocity / gridsize )
    cmax = vpmax * real(DT,kind=CUSTOM_REAL) / distance_min

    ! sets region stability number
    cmax_reg = max(cmax_reg,cmax)

    ! stores for file output
    if (do_output) then
      val_ispec_pmax(ispec) = pmax
      val_ispec_dt(ispec) = dt_max
    endif

  enddo ! ispec
!$OMP ENDDO
!$OMP END PARALLEL

  ! collects for all slices
  vsmin = vsmin_reg
  call min_all_all_cr(vsmin,vsmin_reg)
  vpmax = vpmax_reg
  call max_all_all_cr(vpmax,vpmax_reg)

  dt_max = dt_max_reg
  call min_all_all_cr(dt_max,dt_max_reg)

  pmax = pmax_reg
  call max_all_all_cr(pmax,pmax_reg)

  ! sets global minimum/maximum over all domains
  dt_max_glob = min(dt_max_glob,dt_max_reg)
  pmax_glob = max(pmax_glob,pmax_reg)

  ! region info statistics
  ! main collects info from all MPI slices
  cmax = cmax_reg
  call max_all_cr(cmax,cmax_reg)

  elemsize_min = elemsize_min_reg
  call min_all_cr(elemsize_min,elemsize_min_reg)

  elemsize_max = elemsize_max_reg
  call max_all_cr(elemsize_max,elemsize_max_reg)

  eig_ratio_min = eig_ratio_min_reg
  call min_all_cr(eig_ratio_min,eig_ratio_min_reg)

  eig_ratio_max = eig_ratio_max_reg
  call max_all_cr(eig_ratio_max,eig_ratio_max_reg)

  ! empirical minimum period resolved by mesh
  pmax_empirical = real(T_min_period,kind=CUSTOM_REAL)

  !debug timing
  !tCPU = wtime() - tstart
  !if (myrank == 0) print *,'debug: timing region ',iregion_code,' check mesh resolution = ',sngl(tCPU),'(s)'

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '----------------------------------'
    write(IMAIN,*) '  Verification of mesh parameters:'
    write(IMAIN,*) '----------------------------------'
    write(IMAIN,*) '  Region is ',trim(region(iregion_code))
    write(IMAIN,*)
    if (iregion_code /= IREGION_TRINFINITE .and. iregion_code /= IREGION_INFINITE) then
      if (iregion_code == IREGION_OUTER_CORE) then
        write(IMAIN,*) '  Min Vp = ',vsmin_reg,' (km/s)'
      else
        write(IMAIN,*) '  Min Vs = ',vsmin_reg,' (km/s)'
      endif
      write(IMAIN,*) '  Max Vp = ',vpmax_reg,' (km/s)'
      write(IMAIN,*)
    endif
    write(IMAIN,*) '  Max element edge size = ',elemsize_max_reg,' (km)'
    write(IMAIN,*) '  Min element edge size = ',elemsize_min_reg,' (km)'
    write(IMAIN,*) '  Max/min ratio = ',elemsize_max_reg/elemsize_min_reg
    write(IMAIN,*)
    write(IMAIN,*) '  Max Jacobian eigenvalue ratio = ',eig_ratio_max_reg
    write(IMAIN,*) '  Min Jacobian eigenvalue ratio = ',eig_ratio_min_reg
    write(IMAIN,*)
    if (iregion_code /= IREGION_TRINFINITE .and. iregion_code /= IREGION_INFINITE) then
      write(IMAIN,*) '  Minimum period resolved = ',pmax_reg,' (s)'
      write(IMAIN,*) '  Minimum period resolved (empirical) = ',pmax_empirical,' (s)'
      write(IMAIN,*) '  Maximum suggested time step = ',dt_max_reg,' (s)'
      write(IMAIN,*)
      write(IMAIN,*) '  for DT : ',sngl(DT),' (s)'
      write(IMAIN,*) '  Max stability for wave velocities = ',cmax_reg
    endif
    write(IMAIN,*) '----------------------------------'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! debug
  !if (iregion_code == MAX_NUM_REGIONS) then
  !  ! debug
  !  print *,'*** Maximum suggested time step = ',dt_max_glob
  !endif

  ! debug: saves element flags
  if (do_output) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  saving vtk mesh files for resolution res_minimum_period...'
      call flush_IMAIN()
    endif

    ! minimum period
    filename = prname(1:len_trim(prname))//'res_minimum_period'
    if (USE_VTU_FORMAT) then
      ! binary vtu-format
      call write_VTU_data_elem_cr_binary(nspec,nglob, &
                                         xstore_glob,ystore_glob,zstore_glob, &
                                         ibool,val_ispec_pmax,filename)
    else
      ! vtk-format
      call write_VTK_data_elem_cr(nspec,nglob, &
                                  xstore_glob,ystore_glob,zstore_glob, &
                                  ibool,val_ispec_pmax,filename)
    endif

    ! dt_max
    filename = prname(1:len_trim(prname))//'res_maximum_dt'
    if (USE_VTU_FORMAT) then
      ! binary vtu-format
      call write_VTU_data_elem_cr_binary(nspec,nglob, &
                                         xstore_glob,ystore_glob,zstore_glob, &
                                         ibool,val_ispec_dt,filename)
    else
      ! vtk-format
      call write_VTK_data_elem_cr(nspec,nglob, &
                                  xstore_glob,ystore_glob,zstore_glob, &
                                  ibool,val_ispec_dt,filename)
    endif

    ! free memory
    deallocate(val_ispec_pmax,val_ispec_dt)
  endif

  end subroutine check_mesh_resolution

!
!--------------------------------------------------------------------------------------
!

! note: moved out subroutines of the main routine check_mesh_resolution()
!       to be able to use OpenMP. OpenMP has problems with variables within contained subroutines...


  subroutine get_vpvs_minmax(vpmax,vsmin,ispec,nspec,iregion_code,kappavstore,kappahstore,muvstore,muhstore,rhostore)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ, &
                       PI,GRAV,HUGEVAL,TINYVAL,FOUR_THIRDS, &
                       IREGION_TRINFINITE,IREGION_INFINITE
  use shared_parameters, only: RHOAV,R_PLANET

  !for fully anisotropic models
  !use constants, only: IREGION_CRUST_MANTLE,IREGION_INNER_CORE
  !use meshfem_models_par, only: ANISOTROPIC_INNER_CORE,ANISOTROPIC_3D_MANTLE
  !use regions_mesh_par2, only: &
  !  c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
  !  c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
  !  c36store,c44store,c45store,c46store,c55store,c56store,c66store

  implicit none

  real(kind=CUSTOM_REAL),intent(out) :: vpmax,vsmin
  integer, intent(in) :: ispec,nspec,iregion_code
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: kappavstore,kappahstore, &
    muvstore,muhstore,rhostore

  ! local parameters
  real(kind=CUSTOM_REAL) :: vpv,vph,vsv,vsh
  integer :: i,j,k
  ! scaling factors to re-dimensionalize units
  real(kind=CUSTOM_REAL) :: scaleval

  ! to avoid conversion warnings
  real(kind=CUSTOM_REAL), parameter :: HUGEVAL_REAL = real(HUGEVAL,kind=CUSTOM_REAL)

  ! checks if anything to do for this region
  if (iregion_code == IREGION_TRINFINITE .or. &
      iregion_code == IREGION_INFINITE) then
    ! no material parameters set, return artificial values
    vpmax = 1.0_CUSTOM_REAL
    vsmin = 1.0_CUSTOM_REAL
    ! all done
    return
  endif

  vpmax = - HUGEVAL_REAL
  vsmin = HUGEVAL_REAL

  do k = 1, NGLLZ, NGLLZ-1
    do j = 1, NGLLY, NGLLY-1
      do i = 1, NGLLX, NGLLX-1
        ! vp squared
        ! special cases not required if the crust/mantle routines provide vpv,vph,vsv,vsh already...
        !if (ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
          ! this likely needs to be improved, for now, it just takes the maximum entry of Cij
          !vpv = max(c11store(i,j,k,ispec),c12store(i,j,k,ispec),c13store(i,j,k,ispec),c14store(i,j,k,ispec), &
          !          c15store(i,j,k,ispec),c16store(i,j,k,ispec),c22store(i,j,k,ispec),c23store(i,j,k,ispec), &
          !          c24store(i,j,k,ispec),c25store(i,j,k,ispec),c26store(i,j,k,ispec),c33store(i,j,k,ispec), &
          !          c34store(i,j,k,ispec),c35store(i,j,k,ispec),c36store(i,j,k,ispec),c44store(i,j,k,ispec), &
          !          c45store(i,j,k,ispec),c46store(i,j,k,ispec),c55store(i,j,k,ispec),c56store(i,j,k,ispec), &
          !          c66store(i,j,k,ispec)) / rhostore(i,j,k,ispec)
          ! vph = vpv
          !
          ! we would need to rotate the global cij to local dij by routine rotate_tensor_radial_to_global()
          ! and then use the formula for TI velocities
          ! vph = sqrt(d11/rho)
          ! vpv = sqrt(d33/rho)
          ! vsh = sqrt(d66/rho)
          ! vsv = sqrt(d44/rho)
          ! eta = d13/(d11 - 2.d0*d44)
        !else if (ANISOTROPIC_INNER_CORE .and. iregion_code == IREGION_INNER_CORE) then
          ! this likely needs to be improved, for now, it just takes the maximum entry of Cij (for given symmetry)
          !vpv = max(c11store(i,j,k,ispec),c33store(i,j,k,ispec),c12store(i,j,k,ispec),c13store(i,j,k,ispec), &
          !          c44store(i,j,k,ispec)) / rhostore(i,j,k,ispec)
          !
          ! transversly isotropic
          ! C11 = A = rho * vph**2
          ! C33 = C = rho * vpv**2
          ! C44 = L = rho * vsv**2
          ! C13 = F = eta * (A - 2*L)
          ! C12 = C11 - 2 C66 = A - 2*N = rho * (vph**2 - 2 * vsh**2)
          ! C22 = C11 = A
          ! C23 = C13 = F
          ! C55 = C44 = L
          ! C66 = N = rho * vsh**2 = (C11-C12)/2
          !vpv = c33store(i,j,k,ispec)/rhostore(i,j,k,ispec)  ! vpv squared
          !vph = c11store(i,j,k,ispec)/rhostore(i,j,k,ispec)
        !else
          !vpv = (kappavstore(i,j,k,ispec) + FOUR_THIRDS*muvstore(i,j,k,ispec)) / rhostore(i,j,k,ispec)
          !vph = (kappahstore(i,j,k,ispec) + FOUR_THIRDS*muhstore(i,j,k,ispec)) / rhostore(i,j,k,ispec)
        !endif
        ! here, we assume that the above will be provided by the anisotropic crust/mantle routines and that we have kappav,.. filled
        vpv = (kappavstore(i,j,k,ispec) + FOUR_THIRDS*muvstore(i,j,k,ispec)) / rhostore(i,j,k,ispec)
        vph = (kappahstore(i,j,k,ispec) + FOUR_THIRDS*muhstore(i,j,k,ispec)) / rhostore(i,j,k,ispec)

        ! max
        if (vpv < vph) then
          if (vpmax < vph) vpmax = vph
        else
          if (vpmax < vpv) vpmax = vpv
        endif

        ! vs squared
        ! special cases not required if the crust/mantle routines provide vpv,vph,vsv,vsh already...
        !if (ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
          ! this likely needs to be improved: assumes that vs velocities are specified either by c44,c55 and c66
          !vsv = min(c44store(i,j,k,ispec),c55store(i,j,k,ispec),c66store(i,j,k,ispec)) / rhostore(i,j,k,ispec)
          !vsh = vsv
          !
          ! we would need to rotate the global cij to local dij by routine rotate_tensor_radial_to_global()
          ! and then use the formula for TI velocities
          ! vph = sqrt(d11/rho)
          ! vpv = sqrt(d33/rho)
          ! vsh = sqrt(d66/rho)
          ! vsv = sqrt(d44/rho)
          ! eta = d13/(d11 - 2.d0*d44)
        !else if (ANISOTROPIC_INNER_CORE .and. iregion_code == IREGION_INNER_CORE) then
          ! this likely needs to be improved: assumes that vs velocities are specified by c44 and c66
          !vsv = min(c44store(i,j,k,ispec),c66store(i,j,k,ispec)) / rhostore(i,j,k,ispec)
          ! transversly isotropic
          ! C11 = A = rho * vph**2
          ! C33 = C = rho * vpv**2
          ! C44 = L = rho * vsv**2
          ! C13 = F = eta * (A - 2*L)
          ! C12 = C11 - 2 C66 = A - 2*N = rho * (vph**2 - 2 * vsh**2)
          ! C22 = C11 = A
          ! C23 = C13 = F
          ! C55 = C44 = L
          ! C66 = N = rho * vsh**2 = (C11-C12)/2
          !vsv = c44store(i,j,k,ispec)/rhostore(i,j,k,ispec)  ! vsv squared
          !vsh = 0.5d0*(c11store(i,j,k,ispec)-c12store(i,j,k,ispec))/rhostore(i,j,k,ispec)
        !else
          !vsv = muvstore(i,j,k,ispec) / rhostore(i,j,k,ispec)
          !vsh = muhstore(i,j,k,ispec) / rhostore(i,j,k,ispec)
        !endif
        ! here, we assume that the above will be provided by the anisotropic crust/mantle routines and that we have muv,.. filled
        vsv = muvstore(i,j,k,ispec) / rhostore(i,j,k,ispec)
        vsh = muhstore(i,j,k,ispec) / rhostore(i,j,k,ispec)

        ! ignore fluid regions with Vs = 0 (takes acoustic wave speeds vp)
        if (vsv < TINYVAL) vsv = vpv
        if (vsh < TINYVAL) vsh = vph

        ! min
        if (vsv < vsh) then
          if (vsmin > vsv) vsmin = vsv
        else
          if (vsmin > vsh) vsmin = vsh
        endif
      enddo
    enddo
  enddo

  ! maximum Vp (in km/s)
  scaleval = real(sqrt(PI*GRAV*RHOAV)*(R_PLANET/1000.0d0),kind=CUSTOM_REAL)
  vpmax = sqrt(vpmax) * scaleval
  vsmin = sqrt(vsmin) * scaleval

  end subroutine get_vpvs_minmax

!
!--------------------------------------------------------------------------------------
!

  subroutine get_elem_minmaxsize(elemsize_min,elemsize_max,ispec,nspec,xstore,ystore,zstore)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,HUGEVAL
  use shared_parameters, only: R_PLANET_KM

  implicit none

  real(kind=CUSTOM_REAL),intent(out) :: elemsize_min,elemsize_max
  integer, intent(in) :: ispec,nspec
  double precision,dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: xstore,ystore,zstore

  ! local parameters
  real(kind=CUSTOM_REAL) :: x1,y1,z1,x2,y2,z2,dist
  integer :: i1,i2,j1,j2,k1,k2
  integer :: i,j,k

  ! to avoid conversion warnings
  real(kind=CUSTOM_REAL), parameter :: HUGEVAL_REAL = real(HUGEVAL,kind=CUSTOM_REAL)

  elemsize_min = HUGEVAL_REAL
  elemsize_max = -HUGEVAL_REAL
  ! loops over the four edges that are along X
  i1 = 1
  i2 = NGLLX
  do k = 1, NGLLZ, NGLLZ-1
    do j = 1, NGLLY, NGLLY-1
      x1 = real(xstore(i1,j,k,ispec),kind=CUSTOM_REAL)
      y1 = real(ystore(i1,j,k,ispec),kind=CUSTOM_REAL)
      z1 = real(zstore(i1,j,k,ispec),kind=CUSTOM_REAL)

      x2 = real(xstore(i2,j,k,ispec),kind=CUSTOM_REAL)
      y2 = real(ystore(i2,j,k,ispec),kind=CUSTOM_REAL)
      z2 = real(zstore(i2,j,k,ispec),kind=CUSTOM_REAL)

      dist = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)
      if (elemsize_min > dist) elemsize_min = dist
      if (elemsize_max < dist) elemsize_max = dist
    enddo
  enddo

  ! loops over the four edges that are along Y
  j1 = 1
  j2 = NGLLY
  do k = 1, NGLLZ, NGLLZ-1
    do i = 1, NGLLX, NGLLX-1
      x1 = real(xstore(i,j1,k,ispec),kind=CUSTOM_REAL)
      y1 = real(ystore(i,j1,k,ispec),kind=CUSTOM_REAL)
      z1 = real(zstore(i,j1,k,ispec),kind=CUSTOM_REAL)

      x2 = real(xstore(i,j2,k,ispec),kind=CUSTOM_REAL)
      y2 = real(ystore(i,j2,k,ispec),kind=CUSTOM_REAL)
      z2 = real(zstore(i,j2,k,ispec),kind=CUSTOM_REAL)

      dist = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)
      if (elemsize_min > dist) elemsize_min = dist
      if (elemsize_max < dist) elemsize_max = dist
    enddo
  enddo

  ! loops over the four edges that are along Z
  k1 = 1
  k2 = NGLLZ
  do j = 1, NGLLY, NGLLY-1
    do i = 1, NGLLX, NGLLX-1
      x1 = real(xstore(i,j,k1,ispec),kind=CUSTOM_REAL)
      y1 = real(ystore(i,j,k1,ispec),kind=CUSTOM_REAL)
      z1 = real(zstore(i,j,k1,ispec),kind=CUSTOM_REAL)

      x2 = real(xstore(i,j,k2,ispec),kind=CUSTOM_REAL)
      y2 = real(ystore(i,j,k2,ispec),kind=CUSTOM_REAL)
      z2 = real(zstore(i,j,k2,ispec),kind=CUSTOM_REAL)

      dist = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)
      if (elemsize_min > dist ) elemsize_min = dist
      if (elemsize_max < dist) elemsize_max = dist
    enddo
  enddo
  ! size (in km)
  elemsize_min = sqrt(elemsize_min) * real(R_PLANET_KM,kind=CUSTOM_REAL)
  elemsize_max = sqrt(elemsize_max) * real(R_PLANET_KM,kind=CUSTOM_REAL)

  end subroutine get_elem_minmaxsize

!
!--------------------------------------------------------------------------------------
!


  subroutine get_min_distance_from_second_GLL_points(dx,ispec,nspec,xstore,ystore,zstore)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,HUGEVAL
  use shared_parameters, only: R_PLANET_KM

  implicit none

  real(kind=CUSTOM_REAL),intent(out) :: dx
  integer, intent(in) :: ispec,nspec
  double precision,dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: xstore,ystore,zstore

  ! local parameters
  real(kind=CUSTOM_REAL) :: x1,y1,z1,x2,y2,z2,dist
  integer :: i,j,k

  ! to avoid conversion warnings
  real(kind=CUSTOM_REAL), parameter :: HUGEVAL_REAL = real(HUGEVAL,kind=CUSTOM_REAL)

  dx = HUGEVAL_REAL
  ! loops over the four edges that are along Z
  do k = 1, NGLLZ, NGLLZ-1
    ! front-left
    x1 = real(xstore(1,2,k,ispec),kind=CUSTOM_REAL)
    y1 = real(ystore(1,2,k,ispec),kind=CUSTOM_REAL)
    z1 = real(zstore(1,2,k,ispec),kind=CUSTOM_REAL)

    x2 = real(xstore(2,1,k,ispec),kind=CUSTOM_REAL)
    y2 = real(ystore(2,1,k,ispec),kind=CUSTOM_REAL)
    z2 = real(zstore(2,1,k,ispec),kind=CUSTOM_REAL)

    dist = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)
    if (dx > dist ) dx = dist

    ! rear-left
    x1 = real(xstore(1,NGLLY-1,k,ispec),kind=CUSTOM_REAL)
    y1 = real(ystore(1,NGLLY-1,k,ispec),kind=CUSTOM_REAL)
    z1 = real(zstore(1,NGLLY-1,k,ispec),kind=CUSTOM_REAL)

    x2 = real(xstore(2,NGLLY,k,ispec),kind=CUSTOM_REAL)
    y2 = real(ystore(2,NGLLY,k,ispec),kind=CUSTOM_REAL)
    z2 = real(zstore(2,NGLLY,k,ispec),kind=CUSTOM_REAL)

    dist = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)
    if (dx > dist ) dx = dist

    ! front-right
    x1 = real(xstore(NGLLX,2,k,ispec),kind=CUSTOM_REAL)
    y1 = real(ystore(NGLLX,2,k,ispec),kind=CUSTOM_REAL)
    z1 = real(zstore(NGLLX,2,k,ispec),kind=CUSTOM_REAL)

    x2 = real(xstore(NGLLX-1,1,k,ispec),kind=CUSTOM_REAL)
    y2 = real(ystore(NGLLX-1,1,k,ispec),kind=CUSTOM_REAL)
    z2 = real(zstore(NGLLX-1,1,k,ispec),kind=CUSTOM_REAL)

    dist = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)
    if (dx > dist ) dx = dist

    ! rear-right
    x1 = real(xstore(NGLLX,NGLLY-1,k,ispec),kind=CUSTOM_REAL)
    y1 = real(ystore(NGLLX,NGLLY-1,k,ispec),kind=CUSTOM_REAL)
    z1 = real(zstore(NGLLX,NGLLY-1,k,ispec),kind=CUSTOM_REAL)

    x2 = real(xstore(NGLLX-1,NGLLY,k,ispec),kind=CUSTOM_REAL)
    y2 = real(ystore(NGLLX-1,NGLLY,k,ispec),kind=CUSTOM_REAL)
    z2 = real(zstore(NGLLX-1,NGLLY,k,ispec),kind=CUSTOM_REAL)

    dist = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)
    if (dx > dist ) dx = dist
  enddo

  ! loops over the four edges that are along Y
  do j = 1, NGLLY, NGLLY-1
    ! front-left
    x1 = real(xstore(1,j,2,ispec),kind=CUSTOM_REAL)
    y1 = real(ystore(1,j,2,ispec),kind=CUSTOM_REAL)
    z1 = real(zstore(1,j,2,ispec),kind=CUSTOM_REAL)

    x2 = real(xstore(2,j,1,ispec),kind=CUSTOM_REAL)
    y2 = real(ystore(2,j,1,ispec),kind=CUSTOM_REAL)
    z2 = real(zstore(2,j,1,ispec),kind=CUSTOM_REAL)

    dist = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)
    if (dx > dist ) dx = dist

    ! rear-left
    x1 = real(xstore(1,j,NGLLZ-1,ispec),kind=CUSTOM_REAL)
    y1 = real(ystore(1,j,NGLLZ-1,ispec),kind=CUSTOM_REAL)
    z1 = real(zstore(1,j,NGLLZ-1,ispec),kind=CUSTOM_REAL)

    x2 = real(xstore(2,j,NGLLZ,ispec),kind=CUSTOM_REAL)
    y2 = real(ystore(2,j,NGLLZ,ispec),kind=CUSTOM_REAL)
    z2 = real(zstore(2,j,NGLLZ,ispec),kind=CUSTOM_REAL)

    dist = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)
    if (dx > dist ) dx = dist

    ! front-right
    x1 = real(xstore(NGLLX,j,2,ispec),kind=CUSTOM_REAL)
    y1 = real(ystore(NGLLX,j,2,ispec),kind=CUSTOM_REAL)
    z1 = real(zstore(NGLLX,j,2,ispec),kind=CUSTOM_REAL)

    x2 = real(xstore(NGLLX-1,j,1,ispec),kind=CUSTOM_REAL)
    y2 = real(ystore(NGLLX-1,j,1,ispec),kind=CUSTOM_REAL)
    z2 = real(zstore(NGLLX-1,j,1,ispec),kind=CUSTOM_REAL)

    dist = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)
    if (dx > dist ) dx = dist

    ! rear-right
    x1 = real(xstore(NGLLX,j,NGLLZ-1,ispec),kind=CUSTOM_REAL)
    y1 = real(ystore(NGLLX,j,NGLLZ-1,ispec),kind=CUSTOM_REAL)
    z1 = real(zstore(NGLLX,j,NGLLZ-1,ispec),kind=CUSTOM_REAL)

    x2 = real(xstore(NGLLX-1,j,NGLLZ,ispec),kind=CUSTOM_REAL)
    y2 = real(ystore(NGLLX-1,j,NGLLZ,ispec),kind=CUSTOM_REAL)
    z2 = real(zstore(NGLLX-1,j,NGLLZ,ispec),kind=CUSTOM_REAL)

    dist = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)
    if (dx > dist ) dx = dist
  enddo

  ! loops over the four edges that are along X
  do i = 1, NGLLX, NGLLX-1
    ! front-left
    x1 = real(xstore(i,1,2,ispec),kind=CUSTOM_REAL)
    y1 = real(ystore(i,1,2,ispec),kind=CUSTOM_REAL)
    z1 = real(zstore(i,1,2,ispec),kind=CUSTOM_REAL)

    x2 = real(xstore(i,2,1,ispec),kind=CUSTOM_REAL)
    y2 = real(ystore(i,2,1,ispec),kind=CUSTOM_REAL)
    z2 = real(zstore(i,2,1,ispec),kind=CUSTOM_REAL)

    dist = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)
    if (dx > dist ) dx = dist

    ! rear-left
    x1 = real(xstore(i,1,NGLLZ-1,ispec),kind=CUSTOM_REAL)
    y1 = real(ystore(i,1,NGLLZ-1,ispec),kind=CUSTOM_REAL)
    z1 = real(zstore(i,1,NGLLZ-1,ispec),kind=CUSTOM_REAL)

    x2 = real(xstore(i,2,NGLLZ,ispec),kind=CUSTOM_REAL)
    y2 = real(ystore(i,2,NGLLZ,ispec),kind=CUSTOM_REAL)
    z2 = real(zstore(i,2,NGLLZ,ispec),kind=CUSTOM_REAL)

    dist = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)
    if (dx > dist ) dx = dist

    ! front-right
    x1 = real(xstore(i,NGLLY,2,ispec),kind=CUSTOM_REAL)
    y1 = real(ystore(i,NGLLY,2,ispec),kind=CUSTOM_REAL)
    z1 = real(zstore(i,NGLLY,2,ispec),kind=CUSTOM_REAL)

    x2 = real(xstore(i,NGLLY-1,1,ispec),kind=CUSTOM_REAL)
    y2 = real(ystore(i,NGLLY-1,1,ispec),kind=CUSTOM_REAL)
    z2 = real(zstore(i,NGLLY-1,1,ispec),kind=CUSTOM_REAL)

    dist = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)
    if (dx > dist ) dx = dist

    ! rear-right
    x1 = real(xstore(i,NGLLY,NGLLZ-1,ispec),kind=CUSTOM_REAL)
    y1 = real(ystore(i,NGLLY,NGLLZ-1,ispec),kind=CUSTOM_REAL)
    z1 = real(zstore(i,NGLLY,NGLLZ-1,ispec),kind=CUSTOM_REAL)

    x2 = real(xstore(i,NGLLY-1,NGLLZ,ispec),kind=CUSTOM_REAL)
    y2 = real(ystore(i,NGLLY-1,NGLLZ,ispec),kind=CUSTOM_REAL)
    z2 = real(zstore(i,NGLLY-1,NGLLZ,ispec),kind=CUSTOM_REAL)

    dist = (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2)
    if (dx > dist ) dx = dist
  enddo
  ! size (in km)
  dx = sqrt(dx) * real(R_PLANET_KM,kind=CUSTOM_REAL)

  end subroutine get_min_distance_from_second_GLL_points


!
!--------------------------------------------------------------------------------------
!

  subroutine get_eigenvalues_min_ratio(eig_ratio_min,eig_ratio_max,ispec)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,HUGEVAL

  use regions_mesh_par2, only: &
    xixstore,xiystore,xizstore, &
    etaxstore,etaystore,etazstore, &
    gammaxstore,gammaystore,gammazstore

  implicit none

  real(kind=CUSTOM_REAL),intent(out) :: eig_ratio_min,eig_ratio_max
  integer, intent(in) :: ispec

  ! local parameters
  integer :: i,j,k
  double precision,dimension(3,3) :: jacobian,A
  double precision :: e1,e2,e3
  double precision :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl

  ! to avoid conversion warnings
  real(kind=CUSTOM_REAL), parameter :: HUGEVAL_REAL = real(HUGEVAL,kind=CUSTOM_REAL)

  eig_ratio_min = HUGEVAL_REAL
  eig_ratio_max = - HUGEVAL_REAL

  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        xixl = xixstore(i,j,k,ispec)
        xiyl = xiystore(i,j,k,ispec)
        xizl = xizstore(i,j,k,ispec)
        etaxl = etaxstore(i,j,k,ispec)
        etayl = etaystore(i,j,k,ispec)
        etazl = etazstore(i,j,k,ispec)
        gammaxl = gammaxstore(i,j,k,ispec)
        gammayl = gammaystore(i,j,k,ispec)
        gammazl = gammazstore(i,j,k,ispec)

        ! computes the Jacobian (determinant) (not inversed yet)
        jacobianl = (xixl*(etayl*gammazl-etazl*gammayl) &
                   - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                   + xizl*(etaxl*gammayl-etayl*gammaxl))
        ! checks
        if (jacobianl <= 0.d0 .or. jacobianl /= jacobianl) then
          print *,"Error negative or invalid jacobian in element",ispec,'jacobianl',jacobianl
          stop "Error invalid jacobian in get_eigenvalues_min_ratio()"
        endif

        ! loads the jacobian matrix
        jacobian(1,1) = xixl
        jacobian(2,1) = xiyl
        jacobian(3,1) = xizl
        jacobian(1,2) = etaxl
        jacobian(2,2) = etayl
        jacobian(3,2) = etazl
        jacobian(1,3) = gammaxl
        jacobian(2,3) = gammayl
        jacobian(3,3) = gammazl

        ! since there is a simple algorithm to get eigenvalues of real symmetric matrices,
        ! we build the symmetric matrix: A = J J^t
        !                                then A^t = (J J^t)^t = (J^t)^t (J)^t = J J^t = A
        ! diagonal
        A(1,1) = jacobian(1,1) * jacobian(1,1) + jacobian(1,2) * jacobian(1,2) + jacobian(1,3) * jacobian(1,3)
        A(2,2) = jacobian(2,1) * jacobian(2,1) + jacobian(2,2) * jacobian(2,2) + jacobian(2,3) * jacobian(2,3)
        A(3,3) = jacobian(3,1) * jacobian(3,1) + jacobian(3,2) * jacobian(3,2) + jacobian(3,3) * jacobian(3,3)
        ! upper half
        A(1,2) = jacobian(1,1) * jacobian(2,1) + jacobian(1,2) * jacobian(2,2) + jacobian(1,3) * jacobian(2,3)
        A(1,3) = jacobian(1,1) * jacobian(3,1) + jacobian(1,2) * jacobian(3,2) + jacobian(1,3) * jacobian(3,3)
        A(2,3) = jacobian(2,1) * jacobian(3,1) + jacobian(2,2) * jacobian(3,2) + jacobian(2,3) * jacobian(3,3)

        ! lower half
        !A(2,1) = jacobian(2,1) * jacobian(1,1) + jacobian(2,2) * jacobian(1,2) + jacobian(2,3) * jacobian(1,3)
        A(2,1) = A(1,2)
        A(3,1) = A(1,3)
        A(3,2) = A(2,3)

        ! gets eigenvalues of A (e1 > e2 > e3)
        call get_eigenvalues_sym(A,e1,e2,e3)

        ! smallest/largest eigenvalue of J ( eig(A) = eig(J) * eig(J^t) = eig(J)**2 )
        if (e3 <= 0.d0) stop "Invalid negative smallest eigenvalue in check_mesh_resolution"
        e1 = sqrt(e1)
        e3 = sqrt(e3)

        ! ratio of smallest vs. largest eigenvalue ( == 1 for no distortion)
        if (eig_ratio_min > e3/e1) eig_ratio_min = real(e3/e1,kind=CUSTOM_REAL)
        if (eig_ratio_max < e3/e1) eig_ratio_max = real(e3/e1,kind=CUSTOM_REAL)

      enddo
    enddo
  enddo

  end subroutine get_eigenvalues_min_ratio

!
!--------------------------------------------------------------------------------------
!

  subroutine get_eigenvalues_sym(A,eig1,eig2,eig3)

! modified to Fortran from: https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3.C3.973_matrices
!
! Smith, Oliver K. (April 1961),
! Eigenvalues of a symmetric 3 x 3 matrix,
! communications of the ACM 4 (4): 168, doi:10.1145/355578.366316
!
! Given a real symmetric 3x3 matrix A, compute the eigenvalues
!
! input:
!   A              - symmetric 3x3 matrix
!
! output:
!   eig1,eig2,eig3 - eigenvalues eig1 >= eig2 >= eig3

  implicit none

  double precision,dimension(3,3),intent(in) :: A
  double precision,intent(out) :: eig1,eig2,eig3

  ! local parameters
  double precision :: p1,p2,q,p,phi,r
  double precision :: tr,det
  double precision,dimension(3,3) :: B

  double precision :: detA

  double precision,dimension(3,3), parameter :: Id = reshape((/1,0,0, 0,1,0, 0,0,1/),(/3,3/))
  double precision,parameter :: ONE_THIRD = 1.d0 / 3.d0
  double precision,parameter :: PI = 3.141592653589793d0
  double precision,parameter :: TOL_ZERO = 1.d-30
  ! eigenvalue tolerance (for ocean prem with R_PLANET_KM = 6368, tolerance is ~ 5.3e-5)
  double precision,parameter :: TOL_EIG = 6.d-5

  ! check symmetry
  if (abs(A(1,2) - A(2,1)) > TOL_ZERO .or. abs(A(1,3) - A(3,1)) > TOL_ZERO .or. abs(A(2,3) - A(3,2)) > TOL_ZERO) then
    print *,"Error non-symmetric input matrix A: ",A(:,:)
    stop 'invalid non-symmetric matrix A in check_mesh_resolution'
  endif

  ! det(A)
  detA = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) &
       + A(1,2)*(A(3,1)*A(2,3) - A(2,1)*A(3,3)) &
       + A(1,3)*(A(2,1)*A(3,2) - A(3,1)*A(2,2))

  ! checks
  if (detA <= 0.d0) then
    print *,"Error negative determinant matrix A: ",detA,'matrix',A(:,:)
    stop "invalid input matrix A has negative determinant in check_mesh_resolution"
  endif

  ! power in off-diagonal entries
  p1 = A(1,2)*A(1,2) + A(1,3)*A(1,3) + A(2,3)*A(2,3)
  if (p1 == 0.d0) then
    ! A is diagonal.
    eig1 = A(1,1)
    eig2 = A(2,2)
    eig3 = A(3,3)
    ! sorts eigenvalues
    ! e.g. [ 3 5 8 ] -> [ 5 3 8]
    if (eig1 < eig2) then
      p = eig1
      eig1 = eig2
      eig2 = p
    endif
    ! e.g. [ 5 3 8 ] -> [ 5 8 3 ]
    if (eig2 < eig3) then
      p = eig2
      eig2 = eig3
      eig3 = p
    endif
    ! e.g. [ 5 8 3 ] -> [ 8 5 3 ]
    if (eig1 < eig2) then
      p = eig1
      eig1 = eig2
      eig2 = p
    endif
  else
    tr = A(1,1) + A(2,2) + A(3,3)
    q = ONE_THIRD * tr  ! q = trace(A)/3

    p2 = (A(1,1) - q)**2 + (A(2,2) - q)**2 + (A(3,3) - q)**2 + 2.d0 * p1
    p = sqrt(p2 / 6.d0)

    B(:,:) = (1.d0 / p) * (A(:,:) - q * Id(:,:))       ! Id is the identity matrix

    det = B(1,1)*(B(2,2)*B(3,3) - B(3,2)*B(2,3)) &
        + B(1,2)*(B(3,1)*B(2,3) - B(2,1)*B(3,3)) &
        + B(1,3)*(B(2,1)*B(3,2) - B(3,1)*B(2,2))
    r = det * 0.5d0 ! r = det / 2

    ! In exact arithmetic for a symmetric matrix  -1 <= r <= 1
    ! but computation error can leave it slightly outside this range.
    if (r <= -1.d0) then
      phi = ONE_THIRD * PI
    else if (r >= 1.d0) then
      phi = 0.d0
    else
      phi = ONE_THIRD * acos(r)
    endif

    ! the eigenvalues satisfy eig3 <= eig2 <= eig1
    eig1 = q + 2.d0 * p * cos(phi)
    eig3 = q + 2.d0 * p * cos(phi + 2.d0 * ONE_THIRD * PI)
    eig2 = tr - eig1 - eig3     ! since trace(A) = eig1 + eig2 + eig3

    !print *,'eigenvalue ',eig1,eig2,eig3,'trace / det = ',tr,det
  endif

  ! checks
  if (eig3 == 0.d0) stop 'Invalid zero eigenvalue'
  if (abs(detA - eig1 * eig2 * eig3)/detA > TOL_EIG) then
    print *,'invalid inaccurate eigenvalues:',eig1*eig2*eig3,'instead of ',detA
    stop 'Invalid inaccurate eigenvalues in check_mesh_resolution'
  endif

  end subroutine get_eigenvalues_sym


