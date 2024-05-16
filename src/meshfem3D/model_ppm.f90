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

!--------------------------------------------------------------------------------------------------
!
! PPM - point profile models
!
! for generic models given as depth profiles at lon/lat using a text-file format like:
!
! #lon(deg), lat(deg), depth(km), Vs-perturbation wrt PREM(%), Vs-PREM (km/s)
!  -10.00000       31.00000       40.00000      -1.775005       4.400000
!  -10.00000       32.00000       40.00000      -1.056823       4.400000
! ...
!
!--------------------------------------------------------------------------------------------------

  module model_ppm_par

  ! ----------------------

  ! scale perturbations in shear speed to perturbations in density and vp
  logical,parameter:: SCALE_MODEL = .false.

  ! factor to convert perturbations in shear speed to perturbations in density
  ! taken from s20rts (see also Qin, 2009, sec. 5.2)
  double precision, parameter :: SCALE_RHO = 0.40d0

  ! SCEC version 4 model relationship http://www.data.scec.org/3Dvelocity/
  !double precision, parameter :: SCALE_RHO = 0.254d0

  ! see: P wave seismic velocity and Vp/Vs ratio beneath the Italian peninsula from local earthquake tomography
  ! (Davide Scadi et al.,2008. tectonophysics)
  !! becomes unstable !!
  !double precision, parameter :: SCALE_VP =  1.75d0 !  corresponds to average vp/vs ratio

  ! Zhou et al. 2005: global upper-mantle structure from finite-frequency surface-wave tomography
  ! http://www.gps.caltech.edu/~yingz/pubs/Zhou_JGR_2005.pdf
  !double precision, parameter :: SCALE_VP =  0.5d0 ! by lab measurements Montagner & Anderson, 1989

  ! Qin et al. 2009, sec. 5.2
  double precision, parameter :: SCALE_VP =  0.588d0 ! by Karato, 1993

  ! ----------------------


  ! file
  character(len=*),parameter :: PPM_file_path = "./DATA/PPM/model.txt"

  ! smoothing parameters
  logical,parameter:: GAUSS_SMOOTHING = .false.

  double precision,parameter:: sigma_h = 10.0 ! 50.0  ! km, horizontal
  double precision,parameter:: sigma_v = 10.0 ! 20.0   ! km, vertical

  integer,parameter:: NUM_GAUSSPOINTS = 10

  ! point profile model_variables
  double precision,dimension(:),allocatable :: PPM_dvs,PPM_lat,PPM_lon,PPM_depth

  double precision :: PPM_maxlat,PPM_maxlon,PPM_minlat,PPM_minlon,PPM_maxdepth,PPM_mindepth
  double precision :: PPM_dlat,PPM_dlon,PPM_ddepth,PPM_max_dvs,PPM_min_dvs

  integer :: PPM_num_v,PPM_num_latperlon,PPM_num_lonperdepth

  end module model_ppm_par

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_ppm_broadcast()

! standard routine to setup model

  use constants
  use model_ppm_par

  implicit none

  ! user info
  if (myrank == 0) then
    write(IMAIN,*) 'broadcast model: PPM'
    call flush_IMAIN()
  endif

  ! upper mantle structure
  if (myrank == 0) call read_model_ppm()

  ! broadcast the information read on the main node to all the nodes
  call bcast_all_singlei(PPM_num_v)
  call bcast_all_singlei(PPM_num_latperlon)
  call bcast_all_singlei(PPM_num_lonperdepth)

  if (myrank /= 0) then
    allocate(PPM_lat(PPM_num_v),PPM_lon(PPM_num_v),PPM_depth(PPM_num_v),PPM_dvs(PPM_num_v))
  endif

  call bcast_all_dp(PPM_dvs(1:PPM_num_v),PPM_num_v)
  call bcast_all_dp(PPM_lat(1:PPM_num_v),PPM_num_v)
  call bcast_all_dp(PPM_lon(1:PPM_num_v),PPM_num_v)
  call bcast_all_dp(PPM_depth(1:PPM_num_v),PPM_num_v)

  call bcast_all_singledp(PPM_maxlat)
  call bcast_all_singledp(PPM_minlat)
  call bcast_all_singledp(PPM_maxlon)
  call bcast_all_singledp(PPM_minlon)
  call bcast_all_singledp(PPM_maxdepth)
  call bcast_all_singledp(PPM_mindepth)
  call bcast_all_singledp(PPM_dlat)
  call bcast_all_singledp(PPM_dlon)
  call bcast_all_singledp(PPM_ddepth)

  end subroutine model_ppm_broadcast

!
!--------------------------------------------------------------------------------------------------
!

  subroutine read_model_ppm()

  use constants
  use model_ppm_par

  implicit none

  ! local parameters
  integer ::            ier,counter,i
  double precision ::    lon,lat,depth,dvs,vs
  character(len=256) ::  line

  !e.g. Mediterranean model
  ! counts entries
  counter = 0
  open(unit=IIN,file=trim(PPM_file_path),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    write(IMAIN,*) ' error opening: ',trim(PPM_file_path)
    call flush_IMAIN()
    ! stop
    call exit_mpi(0,"Error opening model ppm")
  endif

  ! first line is text and will be ignored
  read(IIN,'(a256)') line

  ! counts number of data lines
  ier = 0
  do while (ier == 0 )
    read(IIN,*,iostat=ier) lon,lat,depth,dvs,vs
    if (ier == 0) then
      counter = counter + 1
    endif
  enddo
  close(IIN)

  PPM_num_v = counter
  if (counter < 1) then
    write(IMAIN,*)
    write(IMAIN,*) '  model PPM:',PPM_file_path
    write(IMAIN,*) '     no values read in!!!!!!'
    write(IMAIN,*)
    write(IMAIN,*)
    call flush_IMAIN()
    ! stop
    call exit_mpi(0,' no model PPM ')
  else
    write(IMAIN,*)
    write(IMAIN,*) 'model PPM:',trim(PPM_file_path)
    write(IMAIN,*) '  values: ',counter
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  allocate(PPM_lat(counter),PPM_lon(counter),PPM_depth(counter),PPM_dvs(counter))
  PPM_min_dvs = 0.0
  PPM_max_dvs = 0.0
  PPM_dvs(:) = 0.0

  ! vs values
  open(unit=IIN,file=trim(PPM_file_path),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    write(IMAIN,*) ' error opening: ',trim(PPM_file_path)
    call exit_mpi(0,"Error opening model ppm")
  endif
  read(IIN,'(a256)') line   ! first line is text
  counter = 0
  ier = 0
  do while (ier == 0 )
    read(IIN,*,iostat=ier) lon,lat,depth,dvs,vs
    if (ier == 0) then
      counter = counter + 1
      PPM_lat(counter) = lat
      PPM_lon(counter) = lon
      PPM_depth(counter) = depth
      PPM_dvs(counter) = dvs/100.0

      !debug
      !if (abs(depth - 100.0) < 1.e-3) write(IMAIN,*) '  lon/lat/depth : ',lon,lat,depth,' dvs:',dvs
    endif
  enddo
  close(IIN)
  if (counter /= PPM_num_v) then
    write(IMAIN,*)
    write(IMAIN,*) '  model PPM:',PPM_file_path
    write(IMAIN,*) '     error values read in!!!!!!'
    write(IMAIN,*) '  expected: ',PPM_num_v
    write(IMAIN,*) '  got: ',counter
    call exit_mpi(0,' error model PPM ')
  endif

  ! gets depths (in km) of upper and lower limit
  PPM_minlat = minval( PPM_lat(1:PPM_num_v) )
  PPM_maxlat = maxval( PPM_lat(1:PPM_num_v) )

  PPM_minlon = minval( PPM_lon(1:PPM_num_v) )
  PPM_maxlon = maxval( PPM_lon(1:PPM_num_v) )

  PPM_mindepth = minval( PPM_depth(1:PPM_num_v) )
  PPM_maxdepth = maxval( PPM_depth(1:PPM_num_v) )

  PPM_min_dvs = minval(PPM_dvs(1:PPM_num_v))
  PPM_max_dvs = maxval(PPM_dvs(1:PPM_num_v))

  write(IMAIN,*) 'model PPM:'
  write(IMAIN,*) '  latitude min/max   : ',PPM_minlat,PPM_maxlat
  write(IMAIN,*) '  longitude min/max: ',PPM_minlon,PPM_maxlon
  write(IMAIN,*) '  depth min/max      : ',PPM_mindepth,PPM_maxdepth
  write(IMAIN,*)
  write(IMAIN,*) '  dvs min/max : ',PPM_min_dvs,PPM_max_dvs
  write(IMAIN,*)
  if (SCALE_MODEL) then
    write(IMAIN,*) '  scaling: '
    write(IMAIN,*) '    rho: ',SCALE_RHO
    write(IMAIN,*) '    vp : ',SCALE_VP
    write(IMAIN,*)
  endif
  if (GAUSS_SMOOTHING) then
    write(IMAIN,*) '  smoothing: '
    write(IMAIN,*) '    sigma horizontal : ',sigma_h
    write(IMAIN,*) '    sigma vertical   : ',sigma_v
    write(IMAIN,*)
  endif
  call flush_IMAIN()

  ! steps lengths
  PPM_dlat = 0.0d0
  lat = PPM_lat(1)
  do i = 1,PPM_num_v
    if (abs(lat - PPM_lat(i)) > 1.e-15) then
      PPM_dlat = PPM_lat(i) - lat
      exit
    endif
  enddo

  PPM_dlon = 0.0d0
  lon = PPM_lon(1)
  do i = 1,PPM_num_v
    if (abs(lon - PPM_lon(i)) > 1.e-15) then
      PPM_dlon = PPM_lon(i) - lon
      exit
    endif
  enddo

  PPM_ddepth = 0.0d0
  depth = PPM_depth(1)
  do i = 1,PPM_num_v
    if (abs(depth - PPM_depth(i)) > 1.e-15) then
      PPM_ddepth = PPM_depth(i) - depth
      exit
    endif
  enddo

  if (abs(PPM_dlat) < 1.e-15 .or. abs(PPM_dlon) < 1.e-15 .or. abs(PPM_ddepth) < 1.e-15) then
    write(IMAIN,*) '  model PPM:',PPM_file_path
    write(IMAIN,*) '     error in delta values:'
    write(IMAIN,*) '     dlat : ',PPM_dlat,' dlon: ',PPM_dlon,' ddepth: ',PPM_ddepth
    call flush_IMAIN()
    ! stop
    call exit_mpi(0,' error model PPM ')
  else
    write(IMAIN,*) '  model increments:'
    write(IMAIN,*) '  ddepth: ',sngl(PPM_ddepth),' dlat:',sngl(PPM_dlat),' dlon:',sngl(PPM_dlon)
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  PPM_num_latperlon = int( (PPM_maxlat - PPM_minlat) / PPM_dlat) + 1
  PPM_num_lonperdepth = int( (PPM_maxlon - PPM_minlon) / PPM_dlon ) + 1

  end subroutine read_model_ppm


!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_ppm(radius,theta,phi,dvs,dvp,drho)

! returns dvs,dvp and drho for given radius,theta,phi location

  use constants
  use shared_parameters, only: R_PLANET_KM

  use model_ppm_par

  implicit none

  double precision :: radius,theta,phi,dvs,dvp,drho

  ! local parameters
  integer:: i,j,k
  double precision:: lat,lon,r_depth
  double precision:: min_dvs,max_dvs

  double precision:: g_dvs,g_depth,g_lat,g_lon,x,g_weight,weight_sum,weight_prod

  double precision,parameter:: const_a = sigma_v/3.0
  double precision,parameter:: const_b = sigma_h/3.0

  ! initialize
  dvs = 0.0d0
  dvp = 0.0d0
  drho = 0.0d0

  ! depth of given radius (in km)
  r_depth = R_PLANET_KM*(1.0 - radius)  ! radius is normalized between [0,1]
  if (r_depth > PPM_maxdepth .or. r_depth < PPM_mindepth) return

  lat=(PI_OVER_TWO-theta)*RADIANS_TO_DEGREES
  if (lat < PPM_minlat .or. lat > PPM_maxlat ) return

  lon = phi*RADIANS_TO_DEGREES
  if (lon > 180.0d0) lon = lon-360.0d0
  if (lon < PPM_minlon .or. lon > PPM_maxlon ) return

  ! search location value
  if (.not. GAUSS_SMOOTHING) then
    call get_PPMmodel_value(lat,lon,r_depth,dvs)
    return
  endif

  ! loop over neighboring points
  dvs = 0.0
  weight_sum = 0.0

  do i=-NUM_GAUSSPOINTS,NUM_GAUSSPOINTS
    g_depth = r_depth + i*const_a
    do j=-NUM_GAUSSPOINTS,NUM_GAUSSPOINTS
      g_lon = lon + j*const_b/(R_PLANET_KM*DEGREES_TO_RADIANS)
      do k=-NUM_GAUSSPOINTS,NUM_GAUSSPOINTS
        g_lat = lat + k*const_b/(R_PLANET_KM*DEGREES_TO_RADIANS)

        call get_PPMmodel_value(g_lat,g_lon,g_depth,g_dvs)

        ! horizontal weighting
        x = (g_lat-lat)*DEGREES_TO_RADIANS*R_PLANET_KM
        call get_Gaussianweight(x,sigma_h,g_weight)
        g_dvs = g_dvs*g_weight
        weight_prod = g_weight

        x = (g_lon-lon)*DEGREES_TO_RADIANS*R_PLANET_KM
        call get_Gaussianweight(x,sigma_h,g_weight)
        g_dvs = g_dvs*g_weight
        weight_prod = weight_prod * g_weight

        ! vertical weighting
        x = g_depth-r_depth
        call get_Gaussianweight(x,sigma_v,g_weight)
        g_dvs = g_dvs*g_weight
        weight_prod = weight_prod * g_weight

        ! averaging
        weight_sum = weight_sum + weight_prod
        dvs = dvs + g_dvs
      enddo
    enddo
  enddo

  if (weight_sum > 1.e-15) dvs = dvs / weight_sum

  ! store min/max
  max_dvs = PPM_max_dvs
  min_dvs = PPM_min_dvs

  if (dvs > max_dvs ) max_dvs = dvs
  if (dvs < min_dvs ) min_dvs = dvs

  PPM_max_dvs = max_dvs
  PPM_min_dvs = min_dvs

  if (SCALE_MODEL) then
    ! scale density and shear velocity
    drho = SCALE_RHO*dvs
    ! scale vp and shear velocity
    dvp = SCALE_VP*dvs
  endif

  end subroutine model_ppm

!
!--------------------------------------------------------------------------------------------------
!

  subroutine get_PPMmodel_value(lat,lon,depth,dvs)

  use constants
  use model_ppm_par

  implicit none

  double precision lat,lon,depth,dvs

  !integer i,j,k
  !double precision r_top,r_bottom

  integer indexval,num_latperlon,num_lonperdepth

  dvs = 0.0

  if (lat > PPM_maxlat ) return
  if (lat < PPM_minlat ) return
  if (lon > PPM_maxlon ) return
  if (lon < PPM_minlon ) return
  if (depth > PPM_maxdepth ) return
  if (depth < PPM_mindepth ) return

  ! direct access: assumes having a regular interval spacing
  num_latperlon = PPM_num_latperlon ! int( (PPM_maxlat - PPM_minlat) / PPM_dlat) + 1
  num_lonperdepth = PPM_num_lonperdepth ! int( (PPM_maxlon - PPM_minlon) / PPM_dlon ) + 1

  indexval = int( (depth-PPM_mindepth)/PPM_ddepth )*num_lonperdepth*num_latperlon &
           + int( (lon-PPM_minlon)/PPM_dlon )*num_latperlon &
           + int( (lat-PPM_minlat)/PPM_dlat ) + 1
  dvs = PPM_dvs(indexval)

  !  ! loop-wise: slower performance
  !  do i = 1,PPM_num_v
  !    ! depth
  !    r_top = PPM_depth(i)
  !    r_bottom = PPM_depth(i) + PPM_ddepth
  !    if (depth > r_top .and. depth <= r_bottom) then
  !      ! longitude
  !      do j=i,PPM_num_v
  !        if (lon >= PPM_lon(j) .and. lon < PPM_lon(j)+PPM_dlon) then
  !          ! latitude
  !          do k=j,PPM_num_v
  !            if (lat >= PPM_lat(k) .and. lat < PPM_lat(k)+PPM_dlat) then
  !              dvs = PPM_dvs(k)
  !              return
  !            endif
  !          enddo
  !        endif
  !      enddo
  !    endif
  !  enddo

  end subroutine

!
!--------------------------------------------------------------------------------------------------
!

  subroutine get_Gaussianweight(x,sigma,weight)

  use constants

  implicit none

  double precision:: x,sigma,weight

  !double precision,parameter:: one_over2pisqrt = 0.3989422804014327

  ! normalized version
  !weight = one_over2pisqrt*exp(-0.5*x*x/(sigma*sigma))/sigma

  ! only exponential
  weight = exp(-0.5*x*x/(sigma*sigma))

  end subroutine

!
!--------------------------------------------------------------------------------------------------
!

  subroutine smooth_model(nproc_xi,nproc_eta, &
                          rho_vp,rho_vs,nspec_stacey, &
                          iregion_code,xixstore,xiystore,xizstore, &
                          etaxstore,etaystore,etazstore, &
                          gammaxstore,gammaystore,gammazstore, &
                          xstore,ystore,zstore,rhostore, &
                          kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
                          nspec, &
                          NEX_XI,NCHUNKS,ABSORBING_CONDITIONS)

! smooth model parameters

  use constants
  use shared_parameters, only: R_PLANET,R_PLANET_KM

  use model_ppm_par, only: &
    PPM_maxlat,PPM_maxlon,PPM_minlat,PPM_minlon,PPM_maxdepth,PPM_mindepth

  implicit none

  integer, intent(in) :: nproc_xi, nproc_eta

  integer, intent(in) :: NEX_XI

  integer, intent(in) :: nspec,nspec_stacey,NCHUNKS

  logical, intent(in) :: ABSORBING_CONDITIONS

  ! arrays with Jacobian matrix
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(in) :: &
    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore

  ! arrays with mesh parameters
  double precision, intent(in) :: xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision, intent(in) :: ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision, intent(in) :: zstore(NGLLX,NGLLY,NGLLZ,nspec)

  ! for anisotropy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(inout) :: rhostore,kappavstore,kappahstore, &
        muvstore,muhstore,eta_anisostore

  ! Stacey
  real(kind=CUSTOM_REAL), intent(inout) :: rho_vp(NGLLX,NGLLY,NGLLZ,nspec_stacey)
  real(kind=CUSTOM_REAL), intent(inout) :: rho_vs(NGLLX,NGLLY,NGLLZ,nspec_stacey)

  integer, intent(in) :: iregion_code

  ! local parameters
  integer :: i,j,k,ispec

  ! only include the neighboring 3 x 3 slices
  integer, parameter :: NSLICES = 3
  integer ,parameter :: NSLICES2 = NSLICES * NSLICES

  integer :: sizeprocs, ixi, ieta
  integer :: islice(NSLICES2), islice0(NSLICES2), nums

  real(kind=CUSTOM_REAL) :: sigma_h, sigma_h2, sigma_h3, sigma_v, sigma_v2, sigma_v3

  real(kind=CUSTOM_REAL) :: x0, y0, z0, norm, norm_h, norm_v, element_size
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: factor, exp_val

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: jacobian, jacobian0
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xl, yl, zl, xx, yy, zz

  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable :: slice_jacobian
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable :: slice_x, slice_y, slice_z

  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:,:),allocatable :: slice_kernels
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ks_rho,ks_kv,ks_kh,ks_muv,ks_muh,ks_eta,ks_dvp,ks_rhovp,ks_rhovs
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: tk_rho,tk_kv,tk_kh,tk_muv,tk_muh,tk_eta,tk_dvp,tk_rhovp,tk_rhovs
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: bk

  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable:: xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: x, y, z
  real(kind=CUSTOM_REAL), dimension(nspec) :: cx0, cy0, cz0, cx, cy, cz
  double precision :: starttime

  integer :: ii, ispec2, rank, mychunk

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll, wxgll
  double precision, dimension(NGLLY) :: yigll, wygll
  double precision, dimension(NGLLZ) :: zigll, wzgll

  ! array with all the weights in the cube
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: wgll_cube

  real(kind=CUSTOM_REAL), parameter :: ZERO_REAL = 0.0_CUSTOM_REAL

  real(kind=CUSTOM_REAL) :: maxlat,maxlon,maxdepth
  real(kind=CUSTOM_REAL) :: minlat,minlon,mindepth
  real(kind=CUSTOM_REAL) :: radius,theta,phi,lat,lon,r_depth,margin_v,margin_h
  real(kind=CUSTOM_REAL) :: dist_h,dist_v

  double precision,external :: wtime

!----------------------------------------------------------------------------------------------------
  ! smoothing parameters
  logical,parameter :: GAUSS_SMOOTHING = .false. ! set to true to use this smoothing routine

  sigma_h = 100.0_CUSTOM_REAL  ! km, horizontal
  sigma_v = 100.0_CUSTOM_REAL   ! km, vertical

  ! check if smoothing applies
  if (.not. GAUSS_SMOOTHING ) return
!----------------------------------------------------------------------------------------------------

  ! check region: only smooth in mantle & crust
  if (iregion_code /= IREGION_CRUST_MANTLE ) return


  sizeprocs = NCHUNKS*NPROC_XI*NPROC_ETA
  element_size = real((TWO_PI*R_PLANET/1000.d0)/(4*NEX_XI),kind=CUSTOM_REAL)

  if (myrank == 0) then
    write(IMAIN, *) "model smoothing defaults:"
    write(IMAIN, *) "  NPROC_XI , NPROC_ETA, NCHUNKS: ",nproc_xi,nproc_eta,nchunks
    write(IMAIN, *) "  total processors                    : ",sizeprocs
    write(IMAIN, *) "  element size on surface(km): ",element_size
    write(IMAIN, *) "  smoothing sigma horizontal : ",sigma_h," vertical: ", sigma_v
  endif


  if (nchunks == 0) call exit_mpi(myrank,'no chunks')

  element_size = element_size * 1000  ! e.g. 9 km on the surface, 36 km at CMB
  element_size = real(element_size / R_PLANET,kind=CUSTOM_REAL)

  sigma_h = sigma_h * 1000.0_CUSTOM_REAL ! m
  sigma_h = real(sigma_h / R_PLANET,kind=CUSTOM_REAL) ! scale
  sigma_v = sigma_v * 1000.0_CUSTOM_REAL ! m
  sigma_v = real(sigma_v / R_PLANET,kind=CUSTOM_REAL) ! scale

  sigma_h2 = sigma_h ** 2
  sigma_v2 = sigma_v ** 2

  ! search radius
  sigma_h3 = 3.0_CUSTOM_REAL * sigma_h + element_size
  sigma_h3 = sigma_h3 ** 2
  sigma_v3 = 3.0_CUSTOM_REAL * sigma_v + element_size
  sigma_v3 = sigma_v3 ** 2
  ! theoretic normal value
  ! (see integral over -inf to +inf of exp[- x*x/(2*sigma) ] = sigma * sqrt(2*pi) )
  norm_h = real(2.d0*PI*sigma_h**2,kind=CUSTOM_REAL)
  norm_v = real(sqrt(2.d0*PI) * sigma_v,kind=CUSTOM_REAL)
  norm   = norm_h * norm_v

  if (myrank == 0) then
    write(IMAIN, *) "  spectral elements                 : ",nspec
    write(IMAIN, *) "  normalization factor              : ",norm
  endif

  ! GLL points
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)
  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        wgll_cube(i,j,k) = real(wxgll(i)*wygll(j)*wzgll(k),kind=CUSTOM_REAL)
      enddo
    enddo
  enddo

  ! ---- figure out the neighboring 8 or 7 slices: (ichunk,ixi,ieta) index start at 0------
  ! note: ichunk is set to CHUNK_AB etc., while mychunk starts from 0
  mychunk = myrank / (nproc_xi * nproc_eta)
  ieta = (myrank - mychunk * nproc_xi * nproc_eta) / nproc_xi
  ixi = myrank - mychunk * nproc_xi * nproc_eta - ieta * nproc_xi

  ! get the neighboring slices:
  call get_all_eight_slices(mychunk,ixi,ieta, &
        islice0(1),islice0(2),islice0(3),islice0(4),islice0(5),islice0(6),islice0(7),islice0(8), &
        nproc_xi,nproc_eta)

  ! remove the repeated slices (only 8 for corner slices in global case)
  islice(1) = myrank; j = 1
  do i = 1, 8
    if (.not. any(islice(1:i) == islice0(i)) .and. islice0(i) < sizeprocs) then
      j = j + 1
      islice(j) = islice0(i)
    endif
  enddo
  nums = j

  if (myrank == 0) then
    write(IMAIN, *) 'slices:',nums
    write(IMAIN, *) '  ',islice(1:nums)
    write(IMAIN, *)
    call flush_IMAIN()
  endif

  ! read in the topology files of the current and neighboring slices
  ! read in myrank slice
  xl(:,:,:,:) = real(xstore(:,:,:,:),kind=CUSTOM_REAL)
  yl(:,:,:,:) = real(ystore(:,:,:,:),kind=CUSTOM_REAL)
  zl(:,:,:,:) = real(zstore(:,:,:,:),kind=CUSTOM_REAL)

  ! build Jacobian
  allocate(xix(NGLLX,NGLLY,NGLLZ,nspec),xiy(NGLLX,NGLLY,NGLLZ,nspec),xiz(NGLLX,NGLLY,NGLLZ,nspec))
  xix(:,:,:,:) = xixstore(:,:,:,:)
  xiy(:,:,:,:) = xiystore(:,:,:,:)
  xiz(:,:,:,:) = xizstore(:,:,:,:)

  allocate(etax(NGLLX,NGLLY,NGLLZ,nspec),etay(NGLLX,NGLLY,NGLLZ,nspec),etaz(NGLLX,NGLLY,NGLLZ,nspec))
  etax(:,:,:,:) = etaxstore(:,:,:,:)
  etay(:,:,:,:) = etaystore(:,:,:,:)
  etaz(:,:,:,:) = etazstore(:,:,:,:)

  allocate(gammax(NGLLX,NGLLY,NGLLZ,nspec),gammay(NGLLX,NGLLY,NGLLZ,nspec),gammaz(NGLLX,NGLLY,NGLLZ,nspec))
  gammax(:,:,:,:) = gammaxstore(:,:,:,:)
  gammay(:,:,:,:) = gammaystore(:,:,:,:)
  gammaz(:,:,:,:) = gammazstore(:,:,:,:)


  ! get the location of the center of the elements
  do ispec = 1, nspec
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          ! build Jacobian
          !         get derivatives of ux, uy and uz with respect to x, y and z
          xixl = xix(i,j,k,ispec)
          xiyl = xiy(i,j,k,ispec)
          xizl = xiz(i,j,k,ispec)
          etaxl = etax(i,j,k,ispec)
          etayl = etay(i,j,k,ispec)
          etazl = etaz(i,j,k,ispec)
          gammaxl = gammax(i,j,k,ispec)
          gammayl = gammay(i,j,k,ispec)
          gammazl = gammaz(i,j,k,ispec)
          ! compute the Jacobian
          jacobianl = xixl*(etayl*gammazl-etazl*gammayl) - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                        + xizl*(etaxl*gammayl-etayl*gammaxl)

          if (abs(jacobianl) > 1.e-25) then
            jacobianl = 1.0_CUSTOM_REAL / jacobianl
          else
            jacobianl = ZERO_REAL
          endif

          jacobian(i,j,k,ispec) = jacobianl
        enddo
      enddo
    enddo
    cx0(ispec) = (xl(1,1,1,ispec) + xl(NGLLX,NGLLY,NGLLZ,ispec))*0.5
    cy0(ispec) = (yl(1,1,1,ispec) + yl(NGLLX,NGLLY,NGLLZ,ispec))*0.5
    cz0(ispec) = (zl(1,1,1,ispec) + zl(NGLLX,NGLLY,NGLLZ,ispec))*0.5
  enddo
  jacobian0(:,:,:,:) = jacobian(:,:,:,:)

  deallocate(xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

  if (myrank == 0) write(IMAIN, *) 'distributing locations, Jacobians and model values ...'
  call synchronize_all()

  ! get location/Jacobian info from slices
  allocate( slice_x(NGLLX,NGLLY,NGLLZ,NSPEC,nums))
  allocate( slice_y(NGLLX,NGLLY,NGLLZ,NSPEC,nums))
  allocate( slice_z(NGLLX,NGLLY,NGLLZ,NSPEC,nums))
  allocate( slice_jacobian(NGLLX,NGLLY,NGLLZ,NSPEC,nums))
  do rank = 0,sizeprocs-1
    if (rank == myrank) then
      jacobian(:,:,:,:) = jacobian0(:,:,:,:)
      x(:,:,:,:) = real(xstore(:,:,:,:),kind=CUSTOM_REAL)
      y(:,:,:,:) = real(ystore(:,:,:,:),kind=CUSTOM_REAL)
      z(:,:,:,:) = real(zstore(:,:,:,:),kind=CUSTOM_REAL)
    endif
    ! every process broadcasts its info
    call bcast_all_cr(x,NGLLX*NGLLY*NGLLZ*NSPEC)
    call bcast_all_cr(y,NGLLX*NGLLY*NGLLZ*NSPEC)
    call bcast_all_cr(z,NGLLX*NGLLY*NGLLZ*NSPEC)
    call bcast_all_cr(jacobian,NGLLX*NGLLY*NGLLZ*NSPEC)

    ! only relevant process info gets stored
    do ii = 1,nums
      if (islice(ii) == rank) then
        slice_x(:,:,:,:,ii) = x(:,:,:,:)
        slice_y(:,:,:,:,ii) = y(:,:,:,:)
        slice_z(:,:,:,:,ii) = z(:,:,:,:)
        slice_jacobian(:,:,:,:,ii) = jacobian(:,:,:,:)
      endif
    enddo
  enddo

  ! arrays to smooth
  allocate( slice_kernels(NGLLX,NGLLY,NGLLZ,NSPEC,nums,9))
  do rank = 0,sizeprocs-1
    if (rank == myrank) then
      ks_rho(:,:,:,:) = rhostore(:,:,:,:)
      ks_kv(:,:,:,:) = kappavstore(:,:,:,:)
      ks_kh(:,:,:,:) = kappahstore(:,:,:,:)
      ks_muv(:,:,:,:) = muvstore(:,:,:,:)
      ks_muh(:,:,:,:) = muhstore(:,:,:,:)
      ks_eta(:,:,:,:) = eta_anisostore(:,:,:,:)
      if (ABSORBING_CONDITIONS) then
        if (iregion_code == IREGION_CRUST_MANTLE) then
          ks_rhovp(:,:,:,1:nspec_stacey) = rho_vp(:,:,:,1:nspec_stacey)
          ks_rhovs(:,:,:,1:nspec_stacey) = rho_vs(:,:,:,1:nspec_stacey)
        endif
      endif
      ! in case of
      !if (ANISOTROPIC_INNER_CORE .and. iregion_code == IREGION_INNER_CORE) then
      ! or
      !if (ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
      ! or
      !if (ATTENUATION .and. ATTENUATION_3D) then
      ! one should add the c**store and tau_* arrays here as well
    endif
    ! every process broadcasts its info
    call bcast_all_cr(ks_rho,NGLLX*NGLLY*NGLLZ*NSPEC)
    call bcast_all_cr(ks_kv,NGLLX*NGLLY*NGLLZ*NSPEC)
    call bcast_all_cr(ks_kh,NGLLX*NGLLY*NGLLZ*NSPEC)
    call bcast_all_cr(ks_muv,NGLLX*NGLLY*NGLLZ*NSPEC)
    call bcast_all_cr(ks_muh,NGLLX*NGLLY*NGLLZ*NSPEC)
    call bcast_all_cr(ks_eta,NGLLX*NGLLY*NGLLZ*NSPEC)
    call bcast_all_cr(ks_dvp,NGLLX*NGLLY*NGLLZ*NSPEC)
    call bcast_all_cr(ks_rhovp,NGLLX*NGLLY*NGLLZ*NSPEC)
    call bcast_all_cr(ks_rhovs,NGLLX*NGLLY*NGLLZ*NSPEC)

    ! only relevant process info gets stored
    do ii = 1,nums
      if (islice(ii) == rank) then
        slice_kernels(:,:,:,:,ii,1) = ks_rho(:,:,:,:)
        slice_kernels(:,:,:,:,ii,2) = ks_kv(:,:,:,:)
        slice_kernels(:,:,:,:,ii,3) = ks_kh(:,:,:,:)
        slice_kernels(:,:,:,:,ii,4) = ks_muv(:,:,:,:)
        slice_kernels(:,:,:,:,ii,5) = ks_muh(:,:,:,:)
        slice_kernels(:,:,:,:,ii,6) = ks_eta(:,:,:,:)
        slice_kernels(:,:,:,:,ii,7) = ks_dvp(:,:,:,:)
        slice_kernels(:,:,:,:,ii,8) = ks_rhovp(:,:,:,:)
        slice_kernels(:,:,:,:,ii,9) = ks_rhovs(:,:,:,:)
      endif
    enddo
  enddo

  ! get the global maximum value of the original kernel file
  !call synchronize_all()
  !call max_all_cr(maxval(abs(muvstore(:,:,:,:))), max_old)

  if (myrank == 0) write(IMAIN, *) 'start looping over elements and points for smoothing ...'

! loop over all the slices
  tk_rho(:,:,:,:) = 0.0_CUSTOM_REAL
  tk_kh(:,:,:,:) = 0.0_CUSTOM_REAL
  tk_kv(:,:,:,:) = 0.0_CUSTOM_REAL
  tk_muh(:,:,:,:) = 0.0_CUSTOM_REAL
  tk_muv(:,:,:,:) = 0.0_CUSTOM_REAL
  tk_eta(:,:,:,:) = 0.0_CUSTOM_REAL
  tk_dvp(:,:,:,:) = 0.0_CUSTOM_REAL
  tk_rhovp(:,:,:,:) = 0.0_CUSTOM_REAL
  tk_rhovs(:,:,:,:) = 0.0_CUSTOM_REAL

  bk(:,:,:,:) = 0.0_CUSTOM_REAL
  do ii = 1, nums
    if (myrank == 0) starttime = wtime()
    if (myrank == 0) write(IMAIN, *) '  slice number = ', ii

    ! read in the topology, Jacobian, calculate center of elements
    xx(:,:,:,:) = slice_x(:,:,:,:,ii)
    yy(:,:,:,:) = slice_y(:,:,:,:,ii)
    zz(:,:,:,:) = slice_z(:,:,:,:,ii)
    jacobian(:,:,:,:) = slice_jacobian(:,:,:,:,ii)

    ! get the location of the center of the elements
    do ispec2 = 1, nspec
      cx(ispec2) = (xx(1,1,1,ispec2) + xx(NGLLX,NGLLZ,NGLLY,ispec2))*0.5
      cy(ispec2) = (yy(1,1,1,ispec2) + yy(NGLLX,NGLLZ,NGLLY,ispec2))*0.5
      cz(ispec2) = (zz(1,1,1,ispec2) + zz(NGLLX,NGLLZ,NGLLY,ispec2))*0.5
    enddo

    !if (myrank == 0) write(IMAIN, *) '    location:',cx(1),cy(1),cz(1)
    !if (myrank == 0) write(IMAIN, *) '    dist:',(cx(1)-cx0(1))**2+(cy(1)-cy0(1))**2,(cz(1)-cz0(1))**2
    !if (myrank == 0) write(IMAIN, *) '    sigma:',sigma_h3,sigma_v3

    ! array values
    ks_rho(:,:,:,:) = slice_kernels(:,:,:,:,ii,1)
    ks_kv(:,:,:,:) = slice_kernels(:,:,:,:,ii,2)
    ks_kh(:,:,:,:) = slice_kernels(:,:,:,:,ii,3)
    ks_muv(:,:,:,:) = slice_kernels(:,:,:,:,ii,4)
    ks_muh(:,:,:,:) = slice_kernels(:,:,:,:,ii,5)
    ks_eta(:,:,:,:) = slice_kernels(:,:,:,:,ii,6)
    ks_dvp(:,:,:,:) = slice_kernels(:,:,:,:,ii,7)
    ks_rhovp(:,:,:,:) = slice_kernels(:,:,:,:,ii,8)
    ks_rhovs(:,:,:,:) = slice_kernels(:,:,:,:,ii,9)

    ! loop over elements to be smoothed in the current slice
    do ispec = 1, nspec

      if (myrank == 0 .and. mod(ispec,100) == 0 ) write(IMAIN, *) '    ispec ', ispec,' sec:',wtime()-starttime

      ! --- only double loop over the elements in the search radius ---
      do ispec2 = 1, nspec

        ! calculates horizontal and vertical distance between two element centers

        ! vector approximation
        call get_distance_vec(dist_h,dist_v,cx0(ispec),cy0(ispec),cz0(ispec), &
                              cx(ispec2),cy(ispec2),cz(ispec2))

        ! note: distances and sigmah, sigmav are normalized by R_PLANET

        ! checks distance between centers of elements
        if (dist_h > sigma_h3 .or. dist_v > sigma_v3 ) cycle

        ! integration factors
        factor(:,:,:) = jacobian(:,:,:,ispec2) * wgll_cube(:,:,:)

        ! loop over GLL points of the elements in current slice (ispec)
        do k = 1, NGLLZ
          do j = 1, NGLLY
            do i = 1, NGLLX

              ! current point (i,j,k,ispec) location, Cartesian coordinates
              x0 = xl(i,j,k,ispec)
              y0 = yl(i,j,k,ispec)
              z0 = zl(i,j,k,ispec)

              ! calculate weights based on Gaussian smoothing
              call smoothing_weights_vec(x0,y0,z0,sigma_h2,sigma_v2,exp_val, &
                                         xx(:,:,:,ispec2),yy(:,:,:,ispec2),zz(:,:,:,ispec2))

              ! adds GLL integration weights
              exp_val(:,:,:) = exp_val(:,:,:) * factor(:,:,:)

              ! smoothed kernel values
              tk_rho(i,j,k,ispec) = tk_rho(i,j,k,ispec) + sum(exp_val(:,:,:) * ks_rho(:,:,:,ispec2))
              tk_kv(i,j,k,ispec) = tk_kv(i,j,k,ispec) + sum(exp_val(:,:,:) * ks_kv(:,:,:,ispec2))
              tk_kh(i,j,k,ispec) = tk_kh(i,j,k,ispec) + sum(exp_val(:,:,:) * ks_kh(:,:,:,ispec2))
              tk_muv(i,j,k,ispec) = tk_muv(i,j,k,ispec) + sum(exp_val(:,:,:) * ks_muv(:,:,:,ispec2))
              tk_muh(i,j,k,ispec) = tk_muh(i,j,k,ispec) + sum(exp_val(:,:,:) * ks_muh(:,:,:,ispec2))
              tk_eta(i,j,k,ispec) = tk_eta(i,j,k,ispec) + sum(exp_val(:,:,:) * ks_eta(:,:,:,ispec2))
              tk_dvp(i,j,k,ispec) = tk_dvp(i,j,k,ispec) + sum(exp_val(:,:,:) * ks_dvp(:,:,:,ispec2))
              tk_rhovp(i,j,k,ispec) = tk_rhovp(i,j,k,ispec) + sum(exp_val(:,:,:) * ks_rhovp(:,:,:,ispec2))
              tk_rhovs(i,j,k,ispec) = tk_rhovs(i,j,k,ispec) + sum(exp_val(:,:,:) * ks_rhovs(:,:,:,ispec2))

              ! normalization, integrated values of Gaussian smoothing function
              bk(i,j,k,ispec) = bk(i,j,k,ispec) + sum(exp_val(:,:,:))

            enddo
          enddo
        enddo ! (i,j,k)
      enddo ! (ispec2)
    enddo   ! (ispec)
  enddo     ! islice

  if (myrank == 0) write(IMAIN, *) 'Done with integration ...'

  ! gets depths (in km) of upper and lower limit
  maxlat = real(PPM_maxlat,kind=CUSTOM_REAL)
  minlat = real(PPM_minlat,kind=CUSTOM_REAL)

  maxlon = real(PPM_maxlon,kind=CUSTOM_REAL)
  minlon = real(PPM_minlon,kind=CUSTOM_REAL)

  maxdepth = real(PPM_maxdepth,kind=CUSTOM_REAL)
  mindepth = real(PPM_mindepth,kind=CUSTOM_REAL)

  margin_v = real(sigma_v*R_PLANET/1000.d0,kind=CUSTOM_REAL) ! in km
  margin_h = real(sigma_h*R_PLANET/1000.d0 * 180.d0/(R_PLANET_KM*PI),kind=CUSTOM_REAL) ! in degree

  ! computes the smoothed values
  do ispec = 1, nspec

    ! depth of given radius (in km)
    call xyz_2_rthetaphi(cx0(ispec),cy0(ispec),cz0(ispec),radius,theta,phi)
    r_depth = real(R_PLANET_KM - radius*R_PLANET_KM,kind=CUSTOM_REAL)  ! radius is normalized between [0,1]
    if (r_depth >= maxdepth+margin_v .or. r_depth+margin_v < mindepth) cycle

    lat = real((PI/2.0d0-theta)*180.0d0/PI,kind=CUSTOM_REAL)
    if (lat < minlat-margin_h .or. lat > maxlat+margin_h ) cycle

    lon = real(phi*180.0d0/PI,kind=CUSTOM_REAL)
    if (lon > 180.0_CUSTOM_REAL) lon = lon - 360.0_CUSTOM_REAL
    if (lon < minlon-margin_h .or. lon > maxlon+margin_h ) cycle

    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX

          ! check if bk value has an entry
          if (abs(bk(i,j,k,ispec) ) > 1.e-25) then

            ! check if (integrated) normalization value is close to theoretically one
            if (abs(bk(i,j,k,ispec) - norm) > 1.e-3*norm) then ! check the normalization criterion
              print *, 'Problem here --- ', myrank, ispec, i, j, k, bk(i,j,k,ispec), norm
              call exit_mpi(myrank, 'Error computing Gaussian function on the grid')
            endif

            rhostore(i,j,k,ispec) = tk_rho(i,j,k,ispec) / bk(i,j,k,ispec)
            kappavstore(i,j,k,ispec) = tk_kv(i,j,k,ispec) / bk(i,j,k,ispec)
            kappahstore(i,j,k,ispec) = tk_kh(i,j,k,ispec) / bk(i,j,k,ispec)
            muvstore(i,j,k,ispec) = tk_muv(i,j,k,ispec) / bk(i,j,k,ispec)
            muhstore(i,j,k,ispec) = tk_muh(i,j,k,ispec) / bk(i,j,k,ispec)
            eta_anisostore(i,j,k,ispec) = tk_eta(i,j,k,ispec) / bk(i,j,k,ispec)
          endif

        enddo
      enddo
    enddo
  enddo

  if (ABSORBING_CONDITIONS) then
    if (iregion_code == IREGION_CRUST_MANTLE) then
      do ispec = 1, nspec_stacey

        ! depth of given radius (in km)
        call xyz_2_rthetaphi(cx0(ispec),cy0(ispec),cz0(ispec),radius,theta,phi)
        r_depth = real(R_PLANET_KM - radius*R_PLANET_KM,kind=CUSTOM_REAL)  ! radius is normalized between [0,1]
        if (r_depth >= maxdepth+margin_v .or. r_depth+margin_v < mindepth) cycle

        lat = real((PI/2.0d0-theta)*180.0d0/PI,kind=CUSTOM_REAL)
        if (lat < minlat-margin_h .or. lat > maxlat+margin_h ) cycle

        lon = real(phi*180.0d0/PI,kind=CUSTOM_REAL)
        if (lon > 180.0_CUSTOM_REAL) lon = lon - 360.0_CUSTOM_REAL
        if (lon < minlon-margin_h .or. lon > maxlon+margin_h ) cycle

        do k = 1, NGLLZ
          do j = 1, NGLLY
            do i = 1, NGLLX

              ! check if bk value has an entry
              if (abs(bk(i,j,k,ispec) ) > 1.e-25) then
                rho_vp(i,j,k,ispec) = tk_rhovp(i,j,k,ispec)/bk(i,j,k,ispec)
                rho_vs(i,j,k,ispec) = tk_rhovs(i,j,k,ispec)/bk(i,j,k,ispec)
              endif

            enddo
          enddo
        enddo
      enddo
    endif
  endif

  !if (myrank == 0) write(IMAIN, *) 'Maximum data value before smoothing = ', max_old

  ! the maximum value for the smoothed kernel
  !call synchronize_all()
  !call max_all_cr(maxval(abs(muvstore(:,:,:,:))), max_new)

  !if (myrank == 0) then
  !  write(IMAIN, *) 'Maximum data value after smoothing = ', max_new
  !  write(IMAIN, *)
  !endif
  !call synchronize_all()

  end subroutine


