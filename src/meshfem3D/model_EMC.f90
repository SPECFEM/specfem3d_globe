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

! IRIS Earth Model Collaboration (EMC) models
! http://ds.iris.edu/ds/products/emc/
!
! preliminary implementation to support EMC netCDF models that provide vp, vs and density.
! Intended for regional models defining velocity models for crust & upper mantle regions.
!
! version info:
! - April 2023
!   Julien Thurin (original reading routines from https://github.com/thurinj/specfem_nc_tomo/blob/main/read_ncmodel.f90)
!   Daniel Peter (added basic implementation to read/broadcast/interpolate)

module model_emc_par

  use constants, only: CUSTOM_REAL,MAX_STRING_LEN,IMAIN

  use netcdf

  implicit none

  ! velocity model
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: EMC_vs, EMC_vp, EMC_rho

  ! model array dimensions
  integer :: EMC_dims_nx,EMC_dims_ny,EMC_dims_nz

  ! to mask out missing values
  logical, dimension(:,:,:), allocatable :: EMC_mask

  ! grid point coordinates
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: EMC_lat, EMC_lon, EMC_dep
  integer :: EMC_latlen,EMC_lonlen,EMC_deplen

  ! grid extent
  real(kind=CUSTOM_REAL) :: EMC_lat_min,EMC_lat_max
  real(kind=CUSTOM_REAL) :: EMC_lon_min,EMC_lon_max
  real(kind=CUSTOM_REAL) :: EMC_dep_min,EMC_dep_max

  real(kind=CUSTOM_REAL) :: EMC_dlat,EMC_dlon,EMC_ddep  ! grid increments

  logical :: EMC_regular_grid_lat,EMC_regular_grid_lon,EMC_regular_grid_dep
  logical :: EMC_is_regional

  ! units: 1==m, 2==km, 3==m/s, 4==km/s, 5==g/cm^3, 6==kg/cm^3, 7==kg/m^3
  integer :: EMC_dep_unit
  integer :: EMC_vp_unit
  integer :: EMC_vs_unit
  integer :: EMC_rho_unit

  ! positive depth direction 1==up, 2==down
  integer :: EMC_dep_dir

  ! EMC model info
  character(len=MAX_STRING_LEN) :: EMC_model_title, EMC_model_id

  ! enforce simulation mesh to be within EMC model region
  ! If set to .false., we will allow any simulation mesh and (only) embed EMC model values where defined
  ! into the default background model which is set to isotropic PREM.
  logical, parameter :: ENFORCE_EMC_MESH_REGION = .false.

  !------------------------------------------------
  ! Scaling velocities/density
  !------------------------------------------------
  ! note: scaling velocities might be needed if the EMC model misses some parameters.
  !       it is not clear yet, how the scaling should look like, if we want to provide scaling options,
  !       or if partial EMC models are just not supported.
  !
  !       below are some options (similar to the PPM model in file model_ppm.f90)
  !       for now, we don't support scaling yet and require the EMC model provide vp, vs and rho.
  !
  ! scale (perturbations) in shear speed to (perturbations) in density and vp
  logical, parameter:: SCALE_MODEL = .false.

  !! factor to convert (perturbations) in shear speed to (perturbations) in density
  ! taken from s20rts (see also Qin, 2009, sec. 5.2)
  double precision, parameter :: SCALE_RHO = 0.40d0
  ! SCEC version 4 model relationship http://www.data.scec.org/3Dvelocity/
  !double precision, parameter :: SCALE_RHO = 0.254d0

  !! factor to convert (perturbations) in shear speed to (perturbations) in Vp
  ! see: P wave seismic velocity and Vp/Vs ratio beneath the Italian peninsula from local earthquake tomography
  ! (Davide Scadi et al.,2008. tectonophysics)
  !double precision, parameter :: SCALE_VP =  1.75d0 !  corresponds to average vp/vs ratio (becomes unstable!)
  ! Zhou et al. 2005: global upper-mantle structure from finite-frequency surface-wave tomography
  ! http://www.gps.caltech.edu/~yingz/pubs/Zhou_JGR_2005.pdf
  !double precision, parameter :: SCALE_VP =  0.5d0 ! by lab measurements Montagner & Anderson, 1989
  ! Qin et al. 2009, sec. 5.2
  double precision, parameter :: SCALE_VP =  0.588d0 ! by Karato, 1993

  !------------------------------------------------
  ! Parameters for generic reading of NetCDF files
  !------------------------------------------------
  ! Define expected dimension names
  ! IDEA : Later, will be read from Parfile OR looked up in a hardcoded
  ! dictionnary of 'standard' variable names for lat, lon, depths
  ! see mod_dimnames.f90 for reference
  character(len=16), parameter :: latname_exp = 'latitude'
  character(len=16), parameter :: lonname_exp = 'longitude'
  character(len=16), parameter :: depname_exp = 'depth'

  ! Define expected variable names
  ! IDEA : Later, will be read from Parfile
  character(len=16), parameter :: vpname_exp = 'vp'
  character(len=16), parameter :: vsname_exp = 'vs'
  character(len=16), parameter :: rhoname_exp = 'rho'

  character(len=16), parameter :: elevation_exp = 'elevation'

  ! alternate string versions of latitude, longitude and depth for string comparison
  ! Example of how to use this:
  !   if (any(latnames .eq. 'latitude')) then
  !      print *, 'Found latitude'
  character(len=16), dimension(6), parameter :: latnames = (/ character(len=16) :: &
    'latitude','lat','y','ydim','y_dim','y-dim' /)
  character(len=16), dimension(6), parameter :: lonnames = (/ character(len=16) :: &
    'longitude','lon','x','xdim','x_dim','x-dim' /)
  character(len=16), dimension(6), parameter :: depnames = (/ character(len=16) :: &
    'depth','dep','z','zdim','z_dim','z-dim' /)

  ! alternate string versions of vp, vs and rho for string comparison
  ! Alaska model: defines vpfinal,vsfinal,rhofinal for final velocity model
  character(len=16), dimension(6), parameter :: vpnames = (/ character(len=16) :: &
  'vp','vpfinal','VP','vpvar','vp_var','vp-var' /)
  character(len=16), dimension(6), parameter :: vsnames = (/ character(len=16) :: &
  'vs','vsfinal','VS','vsvar','vs_var','vs-var' /)
  character(len=16), dimension(6), parameter :: rhonames = (/ character(len=16) :: &
  'rho','rhofinal','RHO','rhovar','rho_var','rho-var' /)

  ! verbosity
  logical, parameter :: VERBOSE = .false.

contains

! helper tools

  ! --------------------------------
  ! Subroutine to check status of netcdf operation
  ! --------------------------------
  subroutine check_status(status)

  implicit none
  integer, intent (in) :: status

  if (status /= nf90_noerr) then
    print *,'Netcdf error: ',trim(nf90_strerror(status))
    stop 'Netcdf error'
  endif

  end subroutine check_status

  ! --------------------------------
  ! Subroutine to check the global attributes
  ! --------------------------------
  subroutine check_global_attributes(ncid)

  implicit none
  integer, intent(in) :: ncid
  ! local parameters
  integer :: status
  character(len=MAX_STRING_LEN) :: title,id_string
  character(len=100) :: val_string
  real(kind=CUSTOM_REAL) :: val
  real(kind=CUSTOM_REAL) :: lat_min,lat_max,lon_min,lon_max,dep_min,dep_max
  integer :: unit,dir

  if (VERBOSE) print *,'  global attributes: '

  status = nf90_get_att(ncid, NF90_GLOBAL, 'id', id_string)
  if (status == nf90_noerr) then
    if (VERBOSE) print *,'    id   : ',trim(id_string)
  endif

  status = nf90_get_att(ncid, NF90_GLOBAL, 'title', title)
  if (status == nf90_noerr) then
    if (len_trim(title) > 50) title = title(1:50)//'...'
    if (VERBOSE) print *,'    title: ',trim(title)
  endif

  status = nf90_get_att(ncid, NF90_GLOBAL, 'geospatial_lat_min', val)
  if (status == nf90_noerr) then
    lat_min = val
    if (VERBOSE) print *,'    geospatial_lat_min: ',lat_min
  endif

  status = nf90_get_att(ncid, NF90_GLOBAL, 'geospatial_lat_max', val)
  if (status == nf90_noerr) then
    lat_max = val
    if (VERBOSE) print *,'    geospatial_lat_max: ',lat_max
  endif
  !geospatial_lat_units = "degrees_north" ;
  !geospatial_lat_resolution = 0.1f ;

  status = nf90_get_att(ncid, NF90_GLOBAL, 'geospatial_lon_min', val)
  if (status == nf90_noerr) then
    lon_min = val
    if (VERBOSE) print *,'    geospatial_lon_min: ',lon_min
  endif

  status = nf90_get_att(ncid, NF90_GLOBAL, 'geospatial_lon_max', val)
  if (status == nf90_noerr) then
    lon_max = val
    if (VERBOSE) print *,'    geospatial_lon_max: ',lon_max
  endif
  !geospatial_lon_units = "degrees_east" ;
  !geospatial_lon_resolution = 0.2f ;

  status = nf90_get_att(ncid, NF90_GLOBAL, 'geospatial_vertical_min', val)
  if (status == nf90_noerr) then
    dep_min = val
    if (VERBOSE) print *,'    geospatial_vertical_min: ',dep_min
  endif

  status = nf90_get_att(ncid, NF90_GLOBAL, 'geospatial_vertical_max', val)
  if (status == nf90_noerr) then
    dep_max = val
    if (VERBOSE) print *,'    geospatial_vertical_max: ',dep_max
  endif

  status = nf90_get_att(ncid, NF90_GLOBAL, 'geospatial_vertical_units', val_string)
  if (status == nf90_noerr) then
    ! units: 1==m, 2==km, 3==m/s, 4==km/s, 5==g/cm^3, 6==kg/cm^3, 7==kg/m^3
    if (trim(val_string) == 'm') then
      unit = 1
    else if (trim(val_string) == 'km') then
      unit = 2
    else
      unit = 0
    endif
    if (VERBOSE) print *,'    geospatial_vertical_units: ',trim(val_string),' - unit = ',unit
  endif

  status = nf90_get_att(ncid, NF90_GLOBAL, 'geospatial_vertical_positive', val_string)
  if (status == nf90_noerr) then
    if (trim(val_string) == 'up') then
      dir = 1
    else if (trim(val_string) == 'down') then
      dir = 2
    endif
    if (VERBOSE) print *,'    geospatial_vertical_positive: ',trim(val_string),' - dir = ',dir
  endif

  !geospatial_vertical_positive = "down" ;

  if (VERBOSE) print *,''

  ! stores model information
  EMC_model_id = trim(id_string)
  EMC_model_title = trim(title)

  EMC_lat_min = lat_min
  EMC_lat_max = lat_max

  EMC_lon_min = lon_min
  EMC_lon_max = lon_max

  EMC_dep_min = dep_min
  EMC_dep_max = dep_max
  EMC_dep_unit = unit
  EMC_dep_dir = dir

  ! checks depth unit
  ! units: 1==m, 2==km, 3==m/s, 4==km/s, 5==g/cm^3, 6==kg/cm^3, 7==kg/m^3
  if (EMC_dep_unit /= 1 .and. EMC_dep_unit /= 2) then
    print *,'Error: depth unit ',unit,' is invalid, must be either 1==m or 2==km'
    print *,'Please check EMC model file and attribute geospatial_vertical_units ',trim(val_string)
    stop 'EMC depth unit not recognized'
  endif

  end subroutine check_global_attributes


  ! --------------------------------
  ! Subroutine to return dimension ids, names, and lengths
  ! --------------------------------
  subroutine list_dims(ncid, latid, lonid, depid, latlen, lonlen, deplen)

  implicit none
  integer, intent(in) :: ncid
  integer, intent(out) :: latid, lonid, depid
  integer, intent(out) :: latlen, lonlen, deplen
  ! local parameters
  integer :: dimid
  integer :: dimlen
  character(len=100) :: dimname
  integer :: ndims

  ! initializes
  latid = 0; latlen = 0
  lonid = 0; lonlen = 0
  depid = 0; deplen = 0

  if (VERBOSE) print *,'  dimensions: '

  ! Get number of dimensions
  call check_status(nf90_inquire(ncid, ndimensions=ndims))
  if (VERBOSE) print *,'    number of dimensions: ', ndims

  do dimid = 1, ndims
    call check_status(nf90_inquire_dimension(ncid, dimid, dimname, dimlen))
    if (VERBOSE) print *,'    dimension', dimid, ': ', trim(dimname), ' =', dimlen
    ! Assign dimension ids and lengths
    if (trim(dimname) == latname_exp .or. any(latnames == trim(dimname))) then
      latid = dimid
      latlen = dimlen
    elseif (trim(dimname) == lonname_exp .or. any(lonnames == trim(dimname))) then
      lonid = dimid
      lonlen = dimlen
    elseif (trim(dimname) == depname_exp .or. any(depnames == trim(dimname))) then
      depid = dimid
      deplen = dimlen
    endif
  enddo

  if (VERBOSE) print *,''
  if (VERBOSE) print *,'    dim latid = ', latid,' len = ',latlen
  if (VERBOSE) print *,'    dim lonid = ', lonid,' len = ',lonlen
  if (VERBOSE) print *,'    dim depid = ', depid,' len = ',deplen
  if (VERBOSE) print *,''

  end subroutine list_dims

  ! --------------------------------
  ! Subroutine to check that the expected variable names are present
  ! --------------------------------
  subroutine check_varnames(ncid, varid_vp, varid_vs, varid_rho, varid_lat, varid_lon, varid_dep)

  implicit none
  integer, intent(in) :: ncid
  integer, intent(out) :: varid_vp, varid_vs, varid_rho, varid_lat, varid_lon, varid_dep
  ! local parameters
  integer :: varid
  integer :: nvars
  character(len=100) :: varname

  ! initializes
  varid_vp = 0
  varid_vs = 0
  varid_rho = 0

  varid_lat = 0
  varid_lon = 0
  varid_dep = 0

  if (VERBOSE) print *,'  variables: '

  ! Get number of variables
  call check_status(nf90_inquire(ncid, nvariables=nvars))
  if (VERBOSE) print *,'    number of variables: ', nvars

  do varid = 1, nvars
    call check_status(nf90_inquire_variable(ncid, varid, varname))
    if (VERBOSE) print *,'    variable', varid, ': ', trim(varname)
    ! Assign variable ids
    if (trim(varname) == latname_exp .or. any(latnames == trim(varname))) then
      varid_lat = varid
    elseif (trim(varname) == lonname_exp .or. any(lonnames == trim(varname))) then
      varid_lon = varid
    elseif (trim(varname) == depname_exp .or. any(depnames == trim(varname))) then
      varid_dep = varid
    elseif (trim(varname) == vpname_exp .or. any(vpnames == trim(varname))) then
      varid_vp = varid
    elseif (trim(varname) == vsname_exp .or. any(vsnames == trim(varname))) then
      varid_vs = varid
    elseif (trim(varname) == rhoname_exp .or. any(rhonames == trim(varname))) then
      varid_rho = varid
    endif
  enddo

  if (VERBOSE) print *,''
  if (VERBOSE) print *,'    varid_lat  = ', varid_lat
  if (VERBOSE) print *,'    varid_lon  = ', varid_lon
  if (VERBOSE) print *,'    varid_dep  = ', varid_dep
  if (VERBOSE) print *,''
  if (VERBOSE) print *,'    varid_vp   = ', varid_vp
  if (VERBOSE) print *,'    varid_vs   = ', varid_vs
  if (VERBOSE) print *,'    varid_rho  = ', varid_rho
  if (VERBOSE) print *,''

  end subroutine check_varnames

  ! --------------------------------
  ! Subroutine to check the dimension order of a variable
  ! --------------------------------
  subroutine check_dimorder(ncid, varid, latid, lonid, depid, varorderdims)

  implicit none
  integer, intent(in) :: ncid
  integer, intent(in) :: varid
  integer, intent(inout) :: latid, lonid, depid
  integer, intent(inout) :: varorderdims(3)
  ! local parameters
  integer, dimension(:),allocatable :: dimids
  integer :: varndim
  integer :: dimid,index,ier
  character(len=100) :: dimname
  logical :: islat, islon, isdep

  ! initializes
  islat = .false.
  islon = .false.
  isdep = .false.
  varorderdims(:) = 0

  if (VERBOSE) print *,'  vp dimension ordering: '

  ! Get number of dimensions
  call check_status(nf90_inquire_variable(ncid, varid, ndims=varndim))
  if (VERBOSE) print *,'    variable ',varid,': number of dimensions stored in variable = ', varndim

  if (varndim > 0) then
    ! allocates array
    allocate(dimids(varndim),stat=ier)
    if (ier /= 0) stop 'Error allocating dimids array'
    dimids(:) = 0

    ! Get dimension ids
    call check_status(nf90_inquire_variable(ncid, varid, dimids=dimids))

    ! Check that the expected dimensions are present
    index = 0
    do dimid = 1, varndim
      call check_status(nf90_inquire_dimension(ncid, dimids(dimid), dimname))
      if (VERBOSE) print *,'    variable dimension', dimid, ': ', trim(dimname)
      ! Assign dimension ids and lengths
      if (trim(dimname) == latname_exp .or. any(latnames == trim(dimname))) then
        ! varorderdims index
        index = index + 1
        ! checks varorderdims index
        if (index > 3) stop 'Error varoderdims index exceeds 3 dimensions'
        ! Write Variable dimension order to varorderdims
        varorderdims(index) = dimids(dimid)
        islat = .true.
        latid = index
      elseif (trim(dimname) == lonname_exp .or. any(lonnames == trim(dimname))) then
        ! varorderdims index
        index = index + 1
        ! checks varorderdims index
        if (index > 3) stop 'Error varoderdims index exceeds 3 dimensions'
        ! Write Variable dimension order to varorderdims
        varorderdims(index) = dimids(dimid)
        islon = .true.
        lonid = index
      elseif (trim(dimname) == depname_exp .or. any(depnames == trim(dimname))) then
        ! varorderdims index
        index = index + 1
        ! checks varorderdims index
        if (index > 3) stop 'Error varoderdims index exceeds 3 dimensions'
        ! Write Variable dimension order to varorderdims
        varorderdims(index) = dimids(dimid)
        isdep = .true.
        depid = index
      endif
    enddo

    ! free memory
    deallocate(dimids)
  endif

  if (VERBOSE) print *,''
  if (VERBOSE) print *,'    Index order of dimensions for stored values: ', varorderdims
  if (VERBOSE) print *,'    latid ',latid,' var dim order = ',varorderdims(latid)
  if (VERBOSE) print *,'    lonid ',lonid,' var dim order = ',varorderdims(lonid)
  if (VERBOSE) print *,'    depid ',depid,' var dim order = ',varorderdims(depid)
  if (VERBOSE) print *,''

  ! checks if something found
  if (.not. islat .or. .not. islon .or. .not. isdep) then
    print *,'Error: expected dimensions lat/lon/dep are not present in variables'
    stop 'Error netcdf variable dimensions'
  endif

  end subroutine check_dimorder


  ! --------------------------------
  ! Subroutine to check the dimension order of a variable
  ! --------------------------------
  subroutine check_attribute_unit(ncid, varid, unit, direction, missing_val)

  implicit none
  integer, intent(in) :: ncid
  integer, intent(in) :: varid
  integer, intent(inout) :: unit,direction
  real(kind=CUSTOM_REAL), intent(inout) :: missing_val
  ! local parameters
  integer :: num_atts,i
  character(len=100) :: name,value

  ! initializes
  unit = 0
  direction = 0
  missing_val = -1.0

  if (VERBOSE) print *,'  variable units: '

  ! Get number of dimensions
  call check_status(nf90_inquire_variable(ncid, varid, natts=num_atts))
  if (VERBOSE) print *,'    variable ',varid,': number of attributes stored in variable = ', num_atts

  if (num_atts > 0) then
    ! get attribute names
    do i = 1, num_atts
      call check_status(nf90_inq_attname(ncid, varid, i, name))
      if (VERBOSE) print *,'    variable attribute', i, ': ', trim(name)
      ! properties
      if (trim(name) == 'units' .or. trim(name) == '_units') then
        ! assigns unit
        ! units: 1==m, 2==km, 3==m/s, 4==km/s, 5==g/cm^3, 6==kg/cm^3, 7==kg/m^3
        ! get attribute value
        call check_status(nf90_get_att(ncid, varid, name, value))
        if (VERBOSE) print *,'      unit: ',trim(value)
        if (trim(value) == 'm') then
          unit = 1
        else if (trim(value) == 'km') then
          unit = 2
        else if (trim(value) == 'm.s-1') then
          unit = 3
        else if (trim(value) == 'km.s-1') then
          unit = 4
        else if (trim(value) == 'g.cm-3' .or. trim(value) == 'g.cm-1') then
          ! Alaska file defines g.cm-1 which is likely an error as density should be g.cm-3
          unit = 5
        else if (trim(value) == 'kg.cm-3' .or. trim(value) == 'kg.cm-1') then
          unit = 6
        else if (trim(value) == 'kg.m-3' .or. trim(value) == 'kg.m-1') then
          unit = 7
        else
          print *,'Error: unit ',trim(value),' is not recognized yet.'
          print *,'Please ask the developers to help and open a github issue with this error report.'
          stop 'EMC model file with invalid unit'
        endif

      else if (trim(name) == 'positive' .or. trim(name) == '_positive') then
        ! assigns direction
        ! get attribute value
        call check_status(nf90_get_att(ncid, varid, name, value))
        if (VERBOSE) print *,'      positive: ',trim(value)
        if (trim(value) == 'up') then
          ! positive above sealevel, negative depths below sealevel
          direction = 1
        else if (trim(value) == 'down') then
          ! positive towards center of earth
          direction = 2
        endif

      else if (trim(name) == 'missing_value') then
        ! assigns value for invalid entries
        call check_status(nf90_get_att(ncid, varid, name, missing_val))
        if (VERBOSE) print *,'      missing value: ',missing_val
      endif
    enddo
  endif

  if (VERBOSE) print *,''
  if (VERBOSE) print *,'    unit = ',unit,' direction = ',direction,' missing value = ',missing_val
  if (VERBOSE) print *,''

  ! checks if something found
  if (unit == 0) then
    print *,'Error: expected attribute units is not present in variable'
    stop 'Error netcdf variable attribute'
  endif

  end subroutine check_attribute_unit


end module model_emc_par


!
!-------------------------------------------------------------------------------------------
!

  subroutine model_emc_broadcast()

  use constants, only: myrank

  use model_emc_par

  implicit none
  integer :: ier,nx,ny,nz
  integer :: latlen,lonlen,deplen

  ! user info
  if (myrank == 0) then
    write(IMAIN,*) 'broadcast model: EMC model'
    call flush_IMAIN()
  endif

  ! reads in EMC model file
  if (myrank == 0) call read_emc_model()

  call synchronize_all()

  ! broadcasts EMC file info
  call bcast_all_singlecr(EMC_lat_min)
  call bcast_all_singlecr(EMC_lat_max)
  call bcast_all_singlecr(EMC_lon_min)
  call bcast_all_singlecr(EMC_lon_max)
  call bcast_all_singlecr(EMC_dep_min)
  call bcast_all_singlecr(EMC_dep_max)
  call bcast_all_singlel(EMC_is_regional)

  ! broadcasts EMC gridding info
  call bcast_all_singlecr(EMC_dlat)
  call bcast_all_singlecr(EMC_dlon)
  call bcast_all_singlecr(EMC_ddep)

  call bcast_all_singlel(EMC_regular_grid_lat)
  call bcast_all_singlel(EMC_regular_grid_lon)
  call bcast_all_singlel(EMC_regular_grid_dep)

  ! broadcasts grid coordinates
  call bcast_all_singlei(EMC_latlen)
  call bcast_all_singlei(EMC_lonlen)
  call bcast_all_singlei(EMC_deplen)

  latlen = EMC_latlen
  lonlen = EMC_lonlen
  deplen = EMC_deplen

  ! allocate grid arrays for all other processes
  if (myrank /= 0) then
    allocate(EMC_lat(latlen), &
             EMC_lon(lonlen), &
             EMC_dep(deplen),stat=ier)
    if (ier /= 0) stop 'Error allocating lat,lon,dep arrays'
    EMC_lat(:) = 0.0; EMC_lon(:) = 0.0; EMC_dep(:) = 0.0
  endif

  ! broadcast the grid
  call bcast_all_cr(EMC_lat,latlen)
  call bcast_all_cr(EMC_lon,lonlen)
  call bcast_all_cr(EMC_dep,deplen)

  ! broadcasts velocity model values
  call bcast_all_singlei(EMC_dims_nx)
  call bcast_all_singlei(EMC_dims_ny)
  call bcast_all_singlei(EMC_dims_nz)

  nx = EMC_dims_nx
  ny = EMC_dims_ny
  nz = EMC_dims_nz

  ! allocate model arrays for all other processes
  if (myrank /= 0) then
    allocate(EMC_vp(nx,ny,nz), &
             EMC_vs(nx,ny,nz), &
             EMC_rho(nx,ny,nz), stat=ier)
    if (ier /= 0) stop 'Error allocating EMC model arrays'
    EMC_vp(:,:,:) = 0.0; EMC_vs(:,:,:) = 0.0; EMC_rho(:,:,:) = 0.0

    ! mask missing values
    allocate(EMC_mask(nx,ny,nz),stat=ier)
    if (ier /= 0) stop 'Error allocating mask'
    EMC_mask(:,:,:) = .false.
  endif

  ! broadcast the velocity model
  call bcast_all_cr(EMC_vp,nx * ny * nz)
  call bcast_all_cr(EMC_vs,nx * ny * nz)
  call bcast_all_cr(EMC_rho,nx * ny * nz)
  call bcast_all_l(EMC_mask,nx * ny * nz)

  end subroutine model_emc_broadcast

!
!-------------------------------------------------------------------------------------------
!

  subroutine read_emc_model()

  use shared_parameters, only: REGIONAL_MESH_CUTOFF,REGIONAL_MESH_CUTOFF_DEPTH, &
    CENTER_LATITUDE_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH, &
    ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES, &
    NCHUNKS

  use model_emc_par

  implicit none

  ! Define variables
  ! ----------------
  integer :: ncid, status !, varid, dimid ! standard netcdf variables
  !integer :: i, j, k ! looping integers
  integer :: latid, lonid, depid ! dimension ids
  integer :: latlen, lonlen, deplen ! dimension lengths
  integer :: dimlens(3) ! dimension lengths (for array allocation purposes)
  !integer :: nvars, ndims, ngatts, unlimdimid ! netcdf file info
  !integer :: dimids(3) ! dimension ids
  !integer, dimension(:), allocatable :: varids ! variable ids
  integer :: varorderdims(3) ! variable dimension order
  ! Define arrays
  ! -------------
  !real, dimension(:), allocatable :: lat, lon, dep ! lat, lon, dep arrays
  !real, dimension(:,:,:), allocatable :: vp, vs, rho ! model arrays
  integer :: nx,ny,nz

  ! Define strings
  ! --------------
  !character(len=100) :: varname ! variable name
  !character(len=100) :: dimname ! dimension name
  !character(len=100) :: attname ! attribute name
  !character(len=100) :: attvalue ! attribute value
  !character(len=100) :: latname, lonname, depname ! dimension names
  !character(len=100), dimension(:), allocatable :: varnames ! variable names

  ! Define other variables
  ! ----------------------
  !logical :: islat, islon, isdep ! flags for dimension names
  !integer :: dimlen ! dimension length
  !integer :: varndim ! variable number of dimensions
  integer :: varid_vp, varid_vs, varid_rho ! variable ids
  integer :: varid_lat, varid_lon, varid_dep
  integer :: unit_vp,unit_vs,unit_rho,unit_dep
  integer :: dir_dep,dir
  integer :: i,ier
  real(kind=CUSTOM_REAL) :: missing_val_vp,missing_val_vs,missing_val_rho,missing_val_dep
  real(kind=CUSTOM_REAL) :: vp_min,vp_max,vs_min,vs_max,rho_min,rho_max
  real(kind=CUSTOM_REAL) :: dx,dx0
  character (len=MAX_STRING_LEN) :: filename
  ! chunk corners
  double precision :: corners_lat(4),corners_lon(4)
  ! earth radii
  double precision, parameter :: RCMB_ = 3480000.d0
  double precision, parameter :: R_EARTH_ = 6371000.d0

  if (VERBOSE) print *,'reading EMC model:'

  ! Open netcdf file with write access
  ! ----------------
  filename = "./DATA/IRIS_EMC/model.nc"

  ! user output
  write(IMAIN,*) '  file   : ',trim(filename)
  write(IMAIN,*)
  call flush_IMAIN()

  status = nf90_open(filename, NF90_NOWRITE, ncid)

  ! checks file
  if (status /= nf90_noerr) then
    print *,'Error opening EMC model file: ',trim(filename)
    print *,'Please check if file exists...'
    print *
    print *,'Netcdf error: ',trim(nf90_strerror(status))
    stop 'Netcdf error'
  endif

  ! global attributes
  call check_global_attributes(ncid)

  ! converts depth info to km (if needed)
  ! units: 1==m, 2==km, 3==m/s, 4==km/s, 5==g/cm^3, 6==kg/cm^3, 7==kg/m^3
  if (EMC_dep_unit == 1) then
    ! converts to km
    EMC_dep_max = EMC_dep_max / 1000.0
    EMC_dep_min = EMC_dep_min / 1000.0
    EMC_dep_unit = 2                    ! in km
  endif
  ! converts depth direction to be positive down (positive towards center of earth, negative above sealevel)
  if (EMC_dep_dir == 1) then
    ! converts to down
    EMC_dep_max = - EMC_dep_max
    EMC_dep_min = - EMC_dep_min
    EMC_dep_dir = 2                    ! down
  endif

  ! flag for regional/global models
  EMC_is_regional = .true.
  if ((EMC_lat_max - EMC_lat_min == 180.0) .and. &
      (EMC_lon_max - EMC_lon_min == 360.0 )) then
    ! global extent
    EMC_is_regional = .false.
  endif

  ! user output
  write(IMAIN,*) '  id     : ',trim(EMC_model_id)
  write(IMAIN,*) '  title  : ',trim(EMC_model_title)
  write(IMAIN,*)
  write(IMAIN,*) '  range  : latitude  min/max = ',EMC_lat_min,'/',EMC_lat_max
  write(IMAIN,*) '           longitude min/max = ',EMC_lon_min,'/',EMC_lon_max
  write(IMAIN,*) '           depth     min/max = ',EMC_dep_min,'/',EMC_dep_max,'(km)'
  if (EMC_is_regional) then
    write(IMAIN,*) '           regional model'
  else
    write(IMAIN,*) '           global model'
  endif
  write(IMAIN,*)
  call flush_IMAIN()

  ! dimension
  call list_dims(ncid, latid, lonid, depid, latlen, lonlen, deplen)

  ! Check that expected dimensions are present and in correct order
  if (latid == 0) stop 'Error latitude dimension not found'
  if (lonid == 0) stop 'Error longitude dimension not found'
  if (depid == 0) stop 'Error depth dimension not found'

  ! variables
  ! Check that the expected variable names are present
  call check_varnames(ncid, varid_vp, varid_vs, varid_rho, varid_lat, varid_lon, varid_dep)

  ! check
  if (varid_lat == 0) stop 'Error lat array variable not found'
  if (varid_lon == 0) stop 'Error lon array variable not found'
  if (varid_dep == 0) stop 'Error dep array variable not found'
  if (varid_vp == 0) stop 'Error vp array variable not found'
  if (varid_vs == 0) stop 'Error vs array variable not found'
  if (varid_rho == 0) stop 'Error rho array variable not found'

  EMC_latlen = latlen
  EMC_lonlen = lonlen
  EMC_deplen = deplen

  ! allocates model arrays
  allocate(EMC_lat(latlen), &
           EMC_lon(lonlen), &
           EMC_dep(deplen),stat=ier)
  if (ier /= 0) stop 'Error allocating lat,lon,dep arrays'
  EMC_lat(:) = 0.0; EMC_lon(:) = 0.0; EMC_dep(:) = 0.0

  ! Read lat,lon,dep
  call check_status(nf90_get_var(ncid, varid_lat, EMC_lat))
  call check_status(nf90_get_var(ncid, varid_lon, EMC_lon))
  call check_status(nf90_get_var(ncid, varid_dep, EMC_dep))

  ! checks consistency between grid values and header infos
  if (abs(minval(EMC_lat) - EMC_lat_min) > 1.e-15) stop 'Error invalid EMC grid lat min'
  if (abs(maxval(EMC_lat) - EMC_lat_max) > 1.e-15) stop 'Error invalid EMC grid lat max'
  if (abs(minval(EMC_lon) - EMC_lon_min) > 1.e-15) stop 'Error invalid EMC grid lon min'
  if (abs(maxval(EMC_lon) - EMC_lon_max) > 1.e-15) stop 'Error invalid EMC grid lon max'
  if (abs(minval(EMC_dep) - EMC_dep_min) > 1.e-15) stop 'Error invalid EMC grid dep min'
  if (abs(maxval(EMC_dep) - EMC_dep_max) > 1.e-15) stop 'Error invalid EMC grid dep max'

  ! checks up/down direction for depth
  call check_attribute_unit(ncid, varid_dep, unit_dep, dir_dep, missing_val_dep)

  ! converts depth to km
  ! units: 1==m, 2==km, 3==m/s, 4==km/s, 5==g/cm^3, 6==kg/cm^3, 7==kg/m^3
  if (unit_dep == 1) then
    ! converts to km
    EMC_dep(:) = EMC_dep(:) / 1000.d0
    unit_dep = 2    ! in km
  endif
  ! converts depth reference direction to positive being down (positive depth below sealevel, negative depth above)
  if (dir_dep == 1) then
    EMC_dep(:) = - EMC_dep(:)
    dir_dep = 2     ! down
  endif

  ! checks regular gridding
  ! lat & lon are usually regular, depth has often variable grid stepping (finer gridding in crust, coarser in mantle)
  ! lat
  EMC_regular_grid_lat = .true.
  dx0 = EMC_lat(2) - EMC_lat(1)
  do i = 2,latlen-1
    dx = EMC_lat(i+1) - EMC_lat(i)
    if (abs(dx - dx0) > 1.e-4 * dx0) then
      EMC_regular_grid_lat = .false.
      exit
    endif
  enddo
  ! lon
  EMC_regular_grid_lon = .true.
  dx0 = EMC_lon(2) - EMC_lon(1)
  do i = 2,lonlen-1
    dx = EMC_lon(i+1) - EMC_lon(i)
    if (abs(dx - dx0) > 1.e-4 * dx0) then
      EMC_regular_grid_lon = .false.
      exit
    endif
  enddo
  ! depth
  EMC_regular_grid_dep = .true.
  dx0 = EMC_dep(2) - EMC_dep(1)
  do i = 2,deplen-1
    dx = EMC_dep(i+1) - EMC_dep(i)
    if (abs(dx - dx0) > 1.e-4 * dx0) then
      EMC_regular_grid_dep = .false.
      exit
    endif
  enddo

  ! grid spacing
  if (EMC_regular_grid_lat) then
    EMC_dlat = (EMC_lat_max - EMC_lat_min) / (latlen-1)
  else
    EMC_dlat = EMC_lat(2) - EMC_lat(1)
  endif
  if (EMC_regular_grid_lon) then
    EMC_dlon = (EMC_lon_max - EMC_lon_min) / (lonlen-1)
  else
    EMC_dlon = EMC_lon(2) - EMC_lon(1)
  endif
  if (EMC_regular_grid_dep) then
    EMC_ddep = (EMC_dep_max - EMC_dep_min) / (deplen-1)
  else
    EMC_ddep = EMC_dep(2) - EMC_dep(1)
  endif

  ! user output
  write(IMAIN,*) '  grid   : lat  min/max = ', minval(EMC_lat),'/',maxval(EMC_lat)
  write(IMAIN,*) '           lon  min/max = ', minval(EMC_lon),'/',maxval(EMC_lon)
  write(IMAIN,*) '           dep  min/max = ', minval(EMC_dep),'/',maxval(EMC_dep)
  write(IMAIN,*)
  write(IMAIN,*) '           regular grid lat/lon/dep = ',EMC_regular_grid_lat,'/', &
                                                          EMC_regular_grid_lon,'/', &
                                                          EMC_regular_grid_dep
  if (EMC_regular_grid_lat) then
    write(IMAIN,*) '           increments dlat = ',EMC_dlat
  else
    write(IMAIN,*) '           increments dlat = variable (',EMC_dlat,',..)'
  endif
  if (EMC_regular_grid_lon) then
    write(IMAIN,*) '           increments dlon = ',EMC_dlon
  else
    write(IMAIN,*) '           increments dlon = variable (',EMC_dlon,',..)'
  endif
  if (EMC_regular_grid_dep) then
    write(IMAIN,*) '           increments ddep = ',EMC_ddep
  else
    write(IMAIN,*) '           increments ddep = variable (',EMC_ddep,',..)'
  endif
  write(IMAIN,*)
  call flush_IMAIN()

  ! model variables dimensions
  ! Check that the variable's dimensions are in the correct order
  call check_dimorder(ncid, varid_vp, latid, lonid, depid, varorderdims)

  ! set dimension lengths
  dimlens(:) = 0
  dimlens(varorderdims(latid)) = latlen
  dimlens(varorderdims(lonid)) = lonlen
  dimlens(varorderdims(depid)) = deplen

  ! check
  if (any(dimlens == 0)) then
    print *,'Error: vp variable has a zero-length dimension',dimlens(:)
    stop 'Error invalid dimension length'
  endif

  ! Allocate variables vp, vs, rho based on dimension lengths determined above
  ! assumes that vp,vs and rho arrays have all the same dimensions
  nx = dimlens(varorderdims(1))
  ny = dimlens(varorderdims(2))
  nz = dimlens(varorderdims(3))

  ! stores grid dimensions
  EMC_dims_nx = nx
  EMC_dims_ny = ny
  EMC_dims_nz = nz

  ! user output
  write(IMAIN,*) '  size   : number of points lat               = ',latlen
  write(IMAIN,*) '           number of points lon               = ',lonlen
  write(IMAIN,*) '           number of points dep               = ',deplen
  write(IMAIN,*)
  write(IMAIN,*) '           array dimensions nx/ny/nz          = ',nx,'/',ny,'/',nz
  write(IMAIN,*) '           array memory required per process  = ', &
                             sngl(3*dble(nx)*dble(ny)*dble(nz)*dble(CUSTOM_REAL)/1024./1024.),'MB'
  write(IMAIN,*)
  call flush_IMAIN()

  ! allocates model arrays
  allocate(EMC_vp(nx,ny,nz), &
           EMC_vs(nx,ny,nz), &
           EMC_rho(nx,ny,nz),stat=ier)
  if (ier /= 0) stop 'Error allocating vp,vs,rho arrays'
  EMC_vp(:,:,:) = 0.0; EMC_vs(:,:,:) = 0.0; EMC_rho(:,:,:) = 0.0

  ! Set the vill_value of the nc variables (vp, vs, rho) to 0 with
  ! nf90_def_var_fill
  ! NOT WORKING.
  ! call check_status(nf90_def_var_fill(ncid, varid_vp, nf90_fill, 0.0))

  ! Read vp, vs, rho. Replace missing values with 0
  call check_status(nf90_get_var(ncid, varid_vp, EMC_vp))
  call check_status(nf90_get_var(ncid, varid_vs, EMC_vs))
  call check_status(nf90_get_var(ncid, varid_rho, EMC_rho))

  ! gets units and missing values
  call check_attribute_unit(ncid, varid_vp, unit_vp, dir, missing_val_vp)
  call check_attribute_unit(ncid, varid_vs, unit_vs, dir, missing_val_vs)
  call check_attribute_unit(ncid, varid_rho, unit_rho, dir, missing_val_rho)

  ! Close netcdf file
  call check_status(nf90_close(ncid))

  ! mask missing values
  allocate(EMC_mask(nx,ny,nz),stat=ier)
  if (ier /= 0) stop 'Error allocating mask'
  EMC_mask(:,:,:) = .false.

  ! vp
  ! here mask is .false. for missing value, .true. for valid points
  where(EMC_vp == missing_val_vp)
    EMC_mask = .false.
  elsewhere
    EMC_mask = .true.
  end where
  ! min/max without missing values
  vp_min = minval(EMC_vp,mask=EMC_mask)
  vp_max = maxval(EMC_vp,mask=EMC_mask)

  ! vs
  ! here mask is .false. for missing value, .true. for valid points
  where(EMC_vs == missing_val_vs)
    EMC_mask = .false.
  elsewhere
    EMC_mask = .true.
  end where
  ! min/max without missing values
  vs_min = minval(EMC_vs,mask=EMC_mask)
  vs_max = maxval(EMC_vs,mask=EMC_mask)

  ! rho
  ! here mask is .false. for missing value, .true. for valid points
  where(EMC_rho == missing_val_rho)
    EMC_mask = .false.
  elsewhere
    EMC_mask = .true.
  end where
  ! min/max without missing values
  rho_min = minval(EMC_rho,mask=EMC_mask)
  rho_max = maxval(EMC_rho,mask=EMC_mask)

  ! overall mask
  ! opposite value as above: total mask is .true. for missing value to mask out those points
  !                          where at least one of the velocity model values (vp, vs or rho) is missing;
  !                          and mask==.false. for valid points.
  EMC_mask(:,:,:) = .false.
  where(EMC_vp == missing_val_vp) EMC_mask = .true.
  where(EMC_vs == missing_val_vs) EMC_mask = .true.
  where(EMC_rho == missing_val_rho) EMC_mask = .true.

  ! user output
  write(IMAIN,*) '  model  : vp  min/max = ', vp_min,'/',vp_max
  write(IMAIN,*) '           vs  min/max = ', vs_min,'/',vs_max
  write(IMAIN,*) '           rho min/max = ', rho_min,'/',rho_max
  write(IMAIN,*)
  write(IMAIN,*) '           number of missing/incomplete model points is ',count(EMC_mask),' out of ',nx*ny*nz
  write(IMAIN,*) '                                                        ',(100.0*count(EMC_mask))/(nx*ny*nz),'%'
  write(IMAIN,*)
  call flush_IMAIN()

  ! consistency conversions
  ! to have the model parameters always in the same unit (velocities in km/s, density in kg/m^3)

  ! converts density to default kg/m^3
  ! units: 1==m, 2==km, 3==m/s, 4==km/s, 5==g/cm^3, 6==kg/cm^3, 7==kg/m^3
  if (unit_rho == 5) then
    ! converts to kg/m^3
    ! rho [kg/m^3] = rho * 1000 [g/cm^3]
    EMC_rho(:,:,:) = EMC_rho(:,:,:) * 1000.d0
    unit_rho = 7 ! kg/m^3
  else if (unit_rho == 6) then
    ! converts to kg/m^3
    ! rho [kg/m^3] = rho * 1000 [kg/cm^3]
    EMC_rho(:,:,:) = EMC_rho(:,:,:) * 1.d6
    unit_rho = 7 ! kg/m^3
  endif

  ! converts velocity to default m/s
  if (unit_vp == 4) then
    ! converts to m/s
    EMC_vp(:,:,:) = EMC_vp(:,:,:) * 1000.d0
    unit_vp = 3
  endif
  if (unit_vs == 4) then
    ! converts to m/s
    EMC_vs(:,:,:) = EMC_vs(:,:,:) * 1000.d0
    unit_vs = 3
  endif

  ! Par_file region info
  if (NCHUNKS /= 6) then
    ! gets chunck corners lat/lon (lon in range [-180,180])
    call determine_chunk_corners_latlon(CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH, &
                                        ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES, &
                                        corners_lat,corners_lon)
    ! adapts lon range to EMC
    if (EMC_lon_min < -180.0) where(corners_lon(:) > 0.d0) corners_lon = corners_lon - 360.0   ! range [-360,0]
    if (EMC_lon_max > 180.0) where(corners_lon(:) < 0.d0) corners_lon = corners_lon + 360.0    ! range [0,360]

    ! user output
    write(IMAIN,*) '  simulation chunk setup:'
    write(IMAIN,*) '    center (lat,lon) at : (',sngl(CENTER_LATITUDE_IN_DEGREES),sngl(CENTER_LONGITUDE_IN_DEGREES),') geocentric'
    write(IMAIN,*) '    rotation            : ',sngl(GAMMA_ROTATION_AZIMUTH)
    write(IMAIN,*) '    width eta/xi        : ',sngl(ANGULAR_WIDTH_ETA_IN_DEGREES),sngl(ANGULAR_WIDTH_XI_IN_DEGREES)
    write(IMAIN,*)
    write(IMAIN,*) '    corners (lat,lon) at: (',sngl(corners_lat(1)),sngl(corners_lon(1)), &
                                              '),(',sngl(corners_lat(2)),sngl(corners_lon(2)),')'
    write(IMAIN,*) '                          (',sngl(corners_lat(3)),sngl(corners_lon(3)), &
                                              '),(',sngl(corners_lat(4)),sngl(corners_lon(4)),') geographic'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  !! checks compatibility with Par_file settings
  if (ENFORCE_EMC_MESH_REGION) then
    ! regional EMC model
    if (EMC_is_regional) then
      ! global simulations not allowed
      if (NCHUNKS == 6) then
        print *,'Error: Global simulation (NCHUNKS==6) not supported for regional EMC models'
        stop 'EMC regional mesh incompatible with global simulations in Par_file'
      endif

      ! Par_file region must be within EMC grid region
      ! checks corners
      ! lat
      if (EMC_lat_min > minval(corners_lat)) then
        print *,'Error: Par_file defined region too large for EMC model.'
        print *,'Par_file : LATITUDE center = ',CENTER_LATITUDE_IN_DEGREES,' with chunk width = ', ANGULAR_WIDTH_ETA_IN_DEGREES
        print *,'           -> minimum lat = ',CENTER_LATITUDE_IN_DEGREES - 0.5 * ANGULAR_WIDTH_ETA_IN_DEGREES
        print *,'  chunck corners (lat,lon) at: (',corners_lat(1),corners_lon(1),'),(',corners_lat(2),corners_lon(2),')'
        print *,'                               (',corners_lat(3),corners_lon(3),'),(',corners_lat(4),corners_lon(4),')'
        print *,'EMC model: lat min        = ',EMC_lat_min
        print *,'Please choose a smaller region in Par_file for this EMC model.'
        stop 'Invalid region in Par_file for this EMC model.'
      endif
      if (EMC_lat_max < maxval(corners_lat)) then
        print *,'Error: Par_file defined region too large for EMC model.'
        print *,'Par_file : LATITUDE center = ',CENTER_LATITUDE_IN_DEGREES,' with chunk width = ', ANGULAR_WIDTH_ETA_IN_DEGREES
        print *,'           -> maximum lat = ',CENTER_LATITUDE_IN_DEGREES + 0.5 * ANGULAR_WIDTH_ETA_IN_DEGREES
        print *,'  chunck corners (lat lon) at: (',corners_lat(1),corners_lon(1),'),(',corners_lat(2),corners_lon(2),')'
        print *,'                               (',corners_lat(3),corners_lon(3),'),(',corners_lat(4),corners_lon(4),')'
        print *,'EMC model: lat max        = ',EMC_lat_max
        print *,'Please choose a smaller region in Par_file for this EMC model.'
        stop 'Invalid region in Par_file for this EMC model.'
      endif
      ! lon
      if (EMC_lon_min > minval(corners_lon)) then
        print *,'Error: Par_file defined region too large for regional EMC models.'
        print *,'Par_file : LONGITUDE center = ',CENTER_LONGITUDE_IN_DEGREES,' with chunk width = ', ANGULAR_WIDTH_XI_IN_DEGREES
        print *,'           -> minimum lon = ',CENTER_LONGITUDE_IN_DEGREES - 0.5 * ANGULAR_WIDTH_XI_IN_DEGREES
        print *,'  chunck corners (lat,lon) at: (',corners_lat(1),corners_lon(1),'),(',corners_lat(2),corners_lon(2),')'
        print *,'                               (',corners_lat(3),corners_lon(3),'),(',corners_lat(4),corners_lon(4),')'
        print *,'EMC model: lon min        = ',EMC_lon_min
        print *,'Please choose a smaller region in Par_file for this EMC model.'
        stop 'Invalid region in Par_file for this EMC model.'
      endif
      if (EMC_lon_max < maxval(corners_lon)) then
        print *,'Error: Par_file defined region too large for regional EMC models.'
        print *,'Par_file : LONGITUDE center = ',CENTER_LONGITUDE_IN_DEGREES,' with chunk width = ', ANGULAR_WIDTH_XI_IN_DEGREES
        print *,'           -> maximum lon = ',CENTER_LONGITUDE_IN_DEGREES + 0.5 * ANGULAR_WIDTH_XI_IN_DEGREES
        print *,'  chunck corners (lat,lon) at: (',corners_lat(1),corners_lon(1),'),(',corners_lat(2),corners_lon(2),')'
        print *,'                               (',corners_lat(3),corners_lon(3),'),(',corners_lat(4),corners_lon(4),')'
        print *,'EMC model: lon max        = ',EMC_lon_max
        print *,'Please choose a smaller region in Par_file for this EMC model.'
        stop 'Invalid region in Par_file for this EMC model.'
      endif

      ! checks depth
      if (EMC_dep_max < (R_EARTH_ - RCMB_)/1000.0) then
        ! requires shallow mesh
        !todo: we might want to enforce this and mesh according to the EMC model extent, like:
        !        REGIONAL_MESH_CUTOFF = .true.
        !        REGIONAL_MESH_CUTOFF_DEPTH = EMC_dep_max
        !      and broadcast to other processes.
        if (.not. REGIONAL_MESH_CUTOFF) then
          print *,'Error: Par_file for regional EMC model requires a mesh cutoff at depth ',EMC_dep_max
          print *,'Please set REGIONAL_MESH_CUTOFF to .true. and choose an appropriate REGIONAL_MESH_CUTOFF_DEPTH'
          stop 'Invalid regional mesh cutoff in Par_file for EMC regional model'
        endif
        ! checks cutoff depth within EMC depth extent
        if (REGIONAL_MESH_CUTOFF_DEPTH > EMC_dep_max) then
          print *,'Error: regional EMC model requires a mesh cutoff at depth ',EMC_dep_max
          print *,'       current mesh cutoff depth is at ',REGIONAL_MESH_CUTOFF_DEPTH
          print *,'Please set REGIONAL_MESH_CUTOFF to .true. and choose an appropriate REGIONAL_MESH_CUTOFF_DEPTH'
          stop 'Invalid regional mesh cutoff in Par_file for EMC regional model'
        endif
      endif
    endif ! EMC_is_regional
  endif

  end subroutine read_emc_model

!
!-------------------------------------------------------------------------------------------
!

  subroutine model_EMC_crustmantle(iregion_code,r,theta,phi,vpv,vph,vsv,vsh,eta_aniso,rho)

  use constants
  use shared_parameters, only: R_PLANET,R_PLANET_KM,RHOAV

  use model_emc_par

  implicit none

  integer ,intent(in) :: iregion_code

  ! radius     - normalized by globe radius [0,1.x]
  ! theta/phi  - colatitude/longitude in rad (range theta/phi = [0,pi] / [0,2pi]
  double precision, intent(in) :: r,theta,phi

  ! absolute values, not perturbations
  double precision, intent(inout) :: vpv,vph,vsv,vsh,eta_aniso,rho

  ! local parameters
  double precision :: r_depth,lat,lon
  double precision :: vp_iso,vs_iso,rho_iso
  double precision :: vpl,vsl,rhol
  double precision :: scaleval_vel,scaleval_rho
  integer :: index_lat,index_lon,index_dep

  ! positioning/interpolation
  double precision :: interp_val
  double precision :: gamma_interp_x,gamma_interp_y,gamma_interp_z
  double precision :: val1,val2,val3,val4,val5,val6,val7,val8
  integer :: ix,iy,iz,Nx,Ny,Nz

  ! checks if anything to do
  ! so far, EMC models only for crust/mantle are supported
  if (iregion_code /= IREGION_CRUST_MANTLE) return

  ! initializes local parameter values
  vpl = 0.d0
  vsl = 0.d0
  rhol = 0.d0

  ! depth of given radius (in km)
  r_depth = R_PLANET_KM * (1.0 - r)  ! radius is normalized between [0,1]

  ! latitude in degrees in [-90,90]
  lat = (PI_OVER_TWO - theta) * RADIANS_TO_DEGREES
  ! longitude in degrees in [0,360]
  lon = phi * RADIANS_TO_DEGREES

  ! shift lon if EMC longitudes start at negative lon
  if (EMC_lon_min < -180.0) then
    if (lon > 0.0) lon = lon - 360.0d0      ! range [-360,0]
  else if (EMC_lon_min < 0.0) then
    if (lon > 180.0) lon = lon - 360.0d0    ! range [-180,180]
  endif

  ! note: the indexing below would work such that is takes the boundary model point as the associated value.
  !       given the boundary point has valid vp/vs/rho values (and is not masked out because of missing values),
  !       it means that the boundary values are spread out to any simulation region outside of the actual EMC box.
  !
  !       here, the choice is that any point outside of the EMC model box falls back to its background velocity model value,
  !       which by default will be the isotropic PREM.
  !
  !       still, this leads to a problem with the depth in case there is a discrepancy between the elevation of the point
  !       when topography is use and if the EMC model starts at zero depth - which would refer to sea level
  !       or surface topography?

  ! checks if point lies in EMC model region:
  !   we add half of the gridding increments to the limits to avoid numerical artefacts,
  !   mostly for points very close to surface or borders.
  !   indexing below will take care of such points slightly beyond edges.
  if (r_depth < (EMC_dep_min - 0.5*EMC_ddep) .or. r_depth > (EMC_dep_max + 0.5*EMC_ddep)) then
    !if (myrank==0) print*,'debug: lat/lon/dep = ',lat,lon,r_depth,' dep min/max ',EMC_dep_min,EMC_dep_max
    return
  endif
  if (lat < (EMC_lat_min - 0.5*EMC_dlat) .or. lat > (EMC_lat_max + 0.5*EMC_dlat)) then
    !if (myrank==0) print*,'debug: lat/lon/dep = ',lat,lon,r_depth,' lat min/max ',EMC_lat_min,EMC_lat_max
    return
  endif
  if (lon < (EMC_lon_min - 0.5*EMC_dlon) .or. lon > (EMC_lon_max + 0.5*EMC_dlon) ) then
    !if (myrank==0) print*,'debug: lat/lon/dep = ',lat,lon,r_depth,' lon min/max ',EMC_lon_min,EMC_lon_max
    return
  endif

  ! finds grid index
  ! lat
  if (EMC_regular_grid_lat) then
    ! indexing starts at 1,..
    index_lat = floor((lat - EMC_lat_min)/EMC_dlat) + 1
  else
    ! gets index of value in EMC_lat closest to given lat
    index_lat = minloc(abs(EMC_lat(:) - lat),dim=1)
  endif
  ! lon
  if (EMC_regular_grid_lon) then
    ! indexing starts at 1,..
    index_lon = floor((lon - EMC_lon_min)/EMC_dlon) + 1
  else
    ! gets index of value in EMC_lon closest to given lon
    index_lat = minloc(abs(EMC_lon(:) - lon),dim=1)
  endif
  ! depth
  if (EMC_regular_grid_dep) then
    ! indexing starts at 1,..
    index_dep = floor((r_depth - EMC_dep_min)/EMC_ddep) + 1
  else
    ! gets index of value in EMC_dep closest to given dep
    index_dep = minloc(abs(EMC_dep(:) - r_depth),dim=1)
    if (index_dep <= 0) index_dep = 1
    ! determine closest, but smaller value (like floor(..) above) than given depth
    ! for example: r_depth = 30.8 -> EMC_dep == 31 closest for minloc index
    !              we want index with EMC_dep == 30 and use the cell between [30,31]km depth for interpolation
    ! in case depth increases in array: ddep > 0 and index with smaller depth is -1
    ! in case depth decreases in array: ddep < 0 and index with smaller depth is +1
    if (EMC_ddep > 0.0) then
      if (EMC_dep(index_dep) > r_depth) index_dep = index_dep -1
    else if (EMC_ddep < 0.0) then
      if (EMC_dep(index_dep) < r_depth) index_dep = index_dep -1
    endif
  endif

  !debug
  !if (r_depth > 30.0 .and. r_depth < 50.0) &
  !  print*,'debug: lat/lon/dep = ',lat,lon,r_depth,'index lat/lon/dep = ',index_lat,index_lon,index_dep, &
  !         'array value lat/lon/dep',EMC_lat(index_lat),EMC_lon(index_lon),EMC_dep(index_dep)

  !! interpolation
  ! sets array indices
  ! note: arrays EMC_vp(:,:,:),.. might use different index orders for different EMC model files.
  !       For example, Alaska file uses nx==lonlen, ny==latlen, nz==deplen,
  !       but other EMC model files might have a different order.
  !       This gets determines by dimlens() and varorderdims(), here we just check the nx index to decide.
  !       This might need further improvement in future to become more general.
  if (EMC_dims_nx == EMC_lonlen) then
    ! order: lon/lat/dep
    ix = index_lon
    iy = index_lat
    iz = index_dep

    Nx = EMC_lonlen
    Ny = EMC_latlen
    Nz = EMC_deplen

    if (ix >= 1 .and. ix <= Nx) gamma_interp_x = lon - dble(EMC_lon(index_lon))
    if (iy >= 1 .and. iy <= Ny) gamma_interp_y = lat - dble(EMC_lat(index_lat))
    if (iz >= 1 .and. iz <= Nz) gamma_interp_z = r_depth - dble(EMC_dep(index_dep))
  else
    ! order: lat/lon/dep
    ix = index_lat
    iy = index_lon
    iz = index_dep

    Nx = EMC_latlen
    Ny = EMC_lonlen
    Nz = EMC_deplen

    if (ix >= 1 .and. ix <= Nx) gamma_interp_x = lat - dble(EMC_lat(index_lat))
    if (iy >= 1 .and. iy <= Ny) gamma_interp_y = lon - dble(EMC_lon(index_lon))
    if (iz >= 1 .and. iz <= Nz) gamma_interp_z = r_depth - dble(EMC_dep(index_dep))
  endif

  ! checks - this check is not needed, we will consider edge points below
  !if (ix < 0 .or. iy < 0 .or. iz < 0) then
  !  print *,'Error: position has invalid index: '
  !  print *,'  rank                   : ',myrank
  !  print *,'  point index            : ',ix,iy,iz
  !  print *,'  location (lat/lon/dep) : ',sngl(lat),sngl(lon),sngl(r_depth)
  !  print *,'  grid point             : ',EMC_lat(index_lat),EMC_lon(index_lon),EMC_dep(index_dep)
  !  call exit_MPI(myrank,'Error point index in EMC model routine')
  !endif

  ! suppress edge effects for points outside of the model
  if (ix < 1) then
    ix = 1
    gamma_interp_x = 0.d0
  endif
  if (ix > Nx-1) then
    ix = Nx-1
    gamma_interp_x = 1.d0
  endif

  if (iy < 1) then
    iy = 1
    gamma_interp_y = 0.d0
  endif
  if (iy > Ny-1) then
    iy = Ny-1
    gamma_interp_y = 1.d0
  endif

  if (iz < 1) then
    iz = 1
    gamma_interp_z = 0.d0
  endif
  if (iz > Nz-1) then
    iz = Nz-1
    gamma_interp_z = 1.d0
  endif

  if (gamma_interp_x < 0.d0) gamma_interp_x = 0.d0
  if (gamma_interp_x > 1.d0) gamma_interp_x = 1.d0
  if (gamma_interp_y < 0.d0) gamma_interp_y = 0.d0
  if (gamma_interp_y > 1.d0) gamma_interp_y = 1.d0
  if (gamma_interp_z < 0.d0) gamma_interp_z = 0.d0
  if (gamma_interp_z > 1.d0) gamma_interp_z = 1.d0


  ! Voigt average of input (for masked points where model values are missing)
  vp_iso = sqrt( (2.d0*vpv*vpv + vph*vph)/3.d0 )
  vs_iso = sqrt( (2.d0*vsv*vsv + vsh*vsh)/3.d0 )
  rho_iso = rho

  ! input background values are non-dimensionalized already
  ! - here we put dimension back to m/s (and kg/m^3) to mix with EMC model values for interpolation
  vp_iso = vp_iso * R_PLANET * sqrt(PI*GRAV*RHOAV)  ! to m/s
  vs_iso = vs_iso * R_PLANET * sqrt(PI*GRAV*RHOAV)  ! to m/s
  rho_iso = rho_iso * RHOAV                         ! to kg/m3

  ! model vp
  val1 = EMC_vp(ix  ,iy  ,iz  )
  val2 = EMC_vp(ix+1,iy  ,iz  )
  val3 = EMC_vp(ix+1,iy+1,iz  )
  val4 = EMC_vp(ix  ,iy+1,iz  )
  val5 = EMC_vp(ix  ,iy  ,iz+1)
  val6 = EMC_vp(ix+1,iy  ,iz+1)
  val7 = EMC_vp(ix+1,iy+1,iz+1)
  val8 = EMC_vp(ix  ,iy+1,iz+1)

  if (EMC_mask(ix  ,iy  ,iz  )) val1 = vp_iso  ! takes input vp (coming from background model PREM)
  if (EMC_mask(ix+1,iy  ,iz  )) val2 = vp_iso
  if (EMC_mask(ix+1,iy+1,iz  )) val3 = vp_iso
  if (EMC_mask(ix  ,iy+1,iz  )) val4 = vp_iso
  if (EMC_mask(ix  ,iy  ,iz+1)) val5 = vp_iso
  if (EMC_mask(ix+1,iy  ,iz+1)) val6 = vp_iso
  if (EMC_mask(ix+1,iy+1,iz+1)) val7 = vp_iso
  if (EMC_mask(ix  ,iy+1,iz+1)) val8 = vp_iso

  ! interpolation rule
  ! use trilinear interpolation in cell to define perturbation value
  interp_val =  &
       val1 * (1.d0-gamma_interp_x) * (1.d0-gamma_interp_y) * (1.d0-gamma_interp_z) + &
       val2 * gamma_interp_x        * (1.d0-gamma_interp_y) * (1.d0-gamma_interp_z) + &
       val3 * gamma_interp_x        * gamma_interp_y        * (1.d0-gamma_interp_z) + &
       val4 * (1.d0-gamma_interp_x) * gamma_interp_y        * (1.d0-gamma_interp_z) + &
       val5 * (1.d0-gamma_interp_x) * (1.d0-gamma_interp_y) * gamma_interp_z + &
       val6 * gamma_interp_x        * (1.d0-gamma_interp_y) * gamma_interp_z + &
       val7 * gamma_interp_x        * gamma_interp_y        * gamma_interp_z + &
       val8 * (1.d0-gamma_interp_x) * gamma_interp_y        * gamma_interp_z
  ! stores interpolated value
  vpl = interp_val

  ! model vs
  val1 = EMC_vs(ix  ,iy  ,iz  )
  val2 = EMC_vs(ix+1,iy  ,iz  )
  val3 = EMC_vs(ix+1,iy+1,iz  )
  val4 = EMC_vs(ix  ,iy+1,iz  )
  val5 = EMC_vs(ix  ,iy  ,iz+1)
  val6 = EMC_vs(ix+1,iy  ,iz+1)
  val7 = EMC_vs(ix+1,iy+1,iz+1)
  val8 = EMC_vs(ix  ,iy+1,iz+1)

  if (EMC_mask(ix  ,iy  ,iz  )) val1 = vs_iso  ! takes input vs (coming from background model PREM)
  if (EMC_mask(ix+1,iy  ,iz  )) val2 = vs_iso
  if (EMC_mask(ix+1,iy+1,iz  )) val3 = vs_iso
  if (EMC_mask(ix  ,iy+1,iz  )) val4 = vs_iso
  if (EMC_mask(ix  ,iy  ,iz+1)) val5 = vs_iso
  if (EMC_mask(ix+1,iy  ,iz+1)) val6 = vs_iso
  if (EMC_mask(ix+1,iy+1,iz+1)) val7 = vs_iso
  if (EMC_mask(ix  ,iy+1,iz+1)) val8 = vs_iso

  ! interpolation rule
  ! use trilinear interpolation in cell to define perturbation value
  interp_val =  &
       val1 * (1.d0-gamma_interp_x) * (1.d0-gamma_interp_y) * (1.d0-gamma_interp_z) + &
       val2 * gamma_interp_x        * (1.d0-gamma_interp_y) * (1.d0-gamma_interp_z) + &
       val3 * gamma_interp_x        * gamma_interp_y        * (1.d0-gamma_interp_z) + &
       val4 * (1.d0-gamma_interp_x) * gamma_interp_y        * (1.d0-gamma_interp_z) + &
       val5 * (1.d0-gamma_interp_x) * (1.d0-gamma_interp_y) * gamma_interp_z + &
       val6 * gamma_interp_x        * (1.d0-gamma_interp_y) * gamma_interp_z + &
       val7 * gamma_interp_x        * gamma_interp_y        * gamma_interp_z + &
       val8 * (1.d0-gamma_interp_x) * gamma_interp_y        * gamma_interp_z
  ! stores interpolated value
  vsl = interp_val

  ! model rho
  val1 = EMC_rho(ix  ,iy  ,iz  )
  val2 = EMC_rho(ix+1,iy  ,iz  )
  val3 = EMC_rho(ix+1,iy+1,iz  )
  val4 = EMC_rho(ix  ,iy+1,iz  )
  val5 = EMC_rho(ix  ,iy  ,iz+1)
  val6 = EMC_rho(ix+1,iy  ,iz+1)
  val7 = EMC_rho(ix+1,iy+1,iz+1)
  val8 = EMC_rho(ix  ,iy+1,iz+1)

  if (EMC_mask(ix  ,iy  ,iz  )) val1 = rho_iso  ! takes input rho (coming from background model PREM)
  if (EMC_mask(ix+1,iy  ,iz  )) val2 = rho_iso
  if (EMC_mask(ix+1,iy+1,iz  )) val3 = rho_iso
  if (EMC_mask(ix  ,iy+1,iz  )) val4 = rho_iso
  if (EMC_mask(ix  ,iy  ,iz+1)) val5 = rho_iso
  if (EMC_mask(ix+1,iy  ,iz+1)) val6 = rho_iso
  if (EMC_mask(ix+1,iy+1,iz+1)) val7 = rho_iso
  if (EMC_mask(ix  ,iy+1,iz+1)) val8 = rho_iso

  ! interpolation rule
  ! use trilinear interpolation in cell to define perturbation value
  interp_val =  &
       val1 * (1.d0-gamma_interp_x) * (1.d0-gamma_interp_y) * (1.d0-gamma_interp_z) + &
       val2 * gamma_interp_x        * (1.d0-gamma_interp_y) * (1.d0-gamma_interp_z) + &
       val3 * gamma_interp_x        * gamma_interp_y        * (1.d0-gamma_interp_z) + &
       val4 * (1.d0-gamma_interp_x) * gamma_interp_y        * (1.d0-gamma_interp_z) + &
       val5 * (1.d0-gamma_interp_x) * (1.d0-gamma_interp_y) * gamma_interp_z + &
       val6 * gamma_interp_x        * (1.d0-gamma_interp_y) * gamma_interp_z + &
       val7 * gamma_interp_x        * gamma_interp_y        * gamma_interp_z + &
       val8 * (1.d0-gamma_interp_x) * gamma_interp_y        * gamma_interp_z
  ! stores interpolated value
  rhol = interp_val

  !debug
  !if (r_depth > 30.1 .and. r_depth < 50.0) &
  !  print*,'debug: lat/lon/dep = ',lat,lon,r_depth,'vp/vs/rho = ',vpl,vsl,rhol,'iso vp/vs/rho',vp_iso,vs_iso,rho_iso, &
  !          'val ',val1,val2,val3,val4,val5,val6,val7,val8

  ! non-dimensionalize
  scaleval_rho = 1.0d0 / RHOAV                             ! from kg/m3
  scaleval_vel = 1.0d0 / (R_PLANET * sqrt(PI*GRAV*RHOAV))  ! from m/s

  rhol = rhol * scaleval_rho
  vpl = vpl * scaleval_vel
  vsl = vsl * scaleval_vel

  ! returns model values if non-zero
  if (vpl > 1.d-15) then
    ! converts isotropic values to transverse isotropic
    vpv = vpl
    vph = vpl
    vsv = vsl
    vsh = vsl
    eta_aniso = 1.d0
    rho = rhol
  endif

  end subroutine model_EMC_crustmantle
