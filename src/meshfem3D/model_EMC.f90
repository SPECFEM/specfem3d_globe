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
! - Feb 2024
!   Daniel (added interpolation methods to fill missing values, fixed final grid cell interpolation, added tapering)

module model_emc_par

  use constants, only: CUSTOM_REAL,TINYVAL,MAX_STRING_LEN,IMAIN

  use netcdf

  implicit none

  ! velocity model
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: EMC_vs, EMC_vp, EMC_rho

  ! average 1D velocity model
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: EMC_avg_vs, EMC_avg_vp, EMC_avg_rho

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
  integer :: EMC_dep_unit = 0
  integer :: EMC_vp_unit = 0
  integer :: EMC_vs_unit = 0
  integer :: EMC_rho_unit = 0

  ! positive depth direction 1==up, 2==down
  integer :: EMC_dep_dir = 0

  ! depth reference level 1==earth surface, 2==sea level
  integer :: EMC_depth_reference_level = 0

  ! EMC model info
  character(len=MAX_STRING_LEN) :: EMC_model_title, EMC_model_id

  !------------------------------------------------
  ! Model modifications
  !------------------------------------------------
  !
  ! extends the model values beyond its model range
  ! If set to .false., we will put the background velocities to points outside the model range
  ! (default is .true. to extend)
  logical, parameter :: EXTEND_EMC_MODEL_RANGE = .true.

  ! taper off values towards background velocities for points outside the model range
  logical, parameter :: TAPER_OUTSIDE_RANGE = .true.
  ! maximum taper distance (in degrees)
  double precision, parameter :: taper_distance_maximum_hori_in_degrees = 5.d0   ! horizontal taper
  double precision, parameter :: taper_distance_maximum_vert_in_km = 50.d0       ! vertical taper

  ! fill missing values
  !
  ! missing values in data arrays will be replaced by interpolated values between defined neighboring grid points
  ! (default is .true. to fill by interpolated values)
  logical, parameter :: FILL_MISSING_VALUES_WITH_INTERPOLATION = .true.
  ! interpolation method: 1 == Shepard (preferred), 2 == nearest, 3 == bilinear
  integer, parameter :: INTERPOLATION_METHOD = 1
  !
  ! missing values in data arrays will be replaced by the average value at each depth
  ! (default is .false. to use interpolated instead)
  logical, parameter :: FILL_MISSING_VALUES_WITH_AVERAGE = .false.

  ! enforce simulation mesh to be within EMC model region
  !
  ! If set to .false., we will allow any simulation mesh and (only) embed EMC model values where defined
  ! into the default background model which is set to isotropic PREM.
  ! (default is .false. to allow for more flexible, larger meshes)
  logical, parameter :: ENFORCE_EMC_MESH_REGION = .false.

  !------------------------------------------------
  ! Scaling velocities/density
  !------------------------------------------------
  ! note: scaling velocities might be needed if the EMC model misses some parameters.
  !       it is not clear yet, how the scaling should look like, if we want to provide scaling options,
  !       or if partial EMC models are just not supported.
  !
  !       by default, we will use Brocher scaling relations:
  !       Brocher 2005, Empirical Relations between Elastic Wavespeeds and Density in the Earth's Crust, BSSA

  !       note that these scaling relations also would have limitations (range of vs, etc.) in their applicability,
  !       which we ignore for simplicity now.
  !       todo: this can be improved in future version of this routine...
  !
  ! scale missing velocity parameters (Vp-to-Density, Vp-to-Vs)
  logical, parameter :: SCALE_MODEL = .true.

  ! here below are some more simple options (similar to the PPM model in file model_ppm.f90), mentioned for reference.
  !
  !! other factors to convert (perturbations) in shear speed to (perturbations) in density
  ! taken from s20rts
  ! (see also Qin et al. 2009, Reliability of mantle tomography models assessed by spectral element simulation, sec. 5.2)
  !double precision, parameter :: SCALE_RHO = 0.40d0
  ! SCEC CVM-H version model relationship: https://strike.scec.org/scecpedia/CVM-H,
  ! uses Brocher 2005 (Empirical Relations between Elastic Wavespeeds and Density in the Earth's Crust, BSSA) for sediments
  !double precision, parameter :: SCALE_RHO = 0.254d0
  !
  !! other factors to convert (perturbations in) shear speed to (perturbations in) Vp
  ! see: P wave seismic velocity and Vp/Vs ratio beneath the Italian peninsula from local earthquake tomography
  ! (Davide Scadi et al.,2008. tectonophysics)
  !double precision, parameter :: SCALE_VP =  1.75d0 !  corresponds to average vp/vs ratio (becomes unstable!)
  ! Zhou et al. 2005: global upper-mantle structure from finite-frequency surface-wave tomography
  ! http://www.gps.caltech.edu/~yingz/pubs/Zhou_JGR_2005.pdf
  !double precision, parameter :: SCALE_VP =  0.5d0 ! by lab measurements Montagner & Anderson, 1989
  ! (Qin et al. 2009, Reliability of mantle tomography models assessed by spectral element simulation, sec. 5.2)
  !double precision, parameter :: SCALE_VP =  0.588d0 ! by Karato, 1993

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

  ! key words
  ! string versions of latitude, longitude and depth for string comparison
  ! Example of how to use this:
  !   if (any(latnames == 'latitude')) then
  !      print *, 'Found latitude'
  character(len=16), dimension(6), parameter :: latnames = (/ character(len=16) :: &
    'latitude','lat','y','ydim','y_dim','y-dim' /)
  character(len=16), dimension(6), parameter :: lonnames = (/ character(len=16) :: &
    'longitude','lon','x','xdim','x_dim','x-dim' /)
  character(len=16), dimension(6), parameter :: depnames = (/ character(len=16) :: &
    'depth','dep','z','zdim','z_dim','z-dim' /)

  ! string versions of vp, vs and rho for string comparison
  ! Alaska model: defines vpfinal,vsfinal,rhofinal for final velocity model
  character(len=16), dimension(6), parameter :: vpnames = (/ character(len=16) :: &
  'vp','vpfinal','VP','vpvar','vp_var','vp-var' /)
  character(len=16), dimension(6), parameter :: vsnames = (/ character(len=16) :: &
  'vs','vsfinal','VS','vsvar','vs_var','vs-var' /)
  character(len=16), dimension(6), parameter :: rhonames = (/ character(len=16) :: &
  'rho','rhofinal','RHO','rhovar','rho_var','rho-var' /)

  ! string versions of units
  ! length
  character(len=16), dimension(2), parameter :: unitnames_m = (/ character(len=16) :: &
  'm','meter' /)
  character(len=16), dimension(2), parameter :: unitnames_km = (/ character(len=16) :: &
  'km','kilometer' /)
  ! density
  character(len=16), dimension(2), parameter :: unitnames_gcm = (/ character(len=16) :: &
  'g.cm-3','g.cm-1' /)  ! Alaska file defines g.cm-1 which is likely an error as density should be g.cm-3
  character(len=16), dimension(2), parameter :: unitnames_kgcm = (/ character(len=16) :: &
  'kg.cm-3','kg.cm-1' /)
  character(len=16), dimension(2), parameter :: unitnames_kgm = (/ character(len=16) :: &
  'kg.m-3','kg.m-1' /)
  ! velocity
  character(len=16), dimension(4), parameter :: unitnames_velms = (/ character(len=16) :: &
  'm.s-1','m/s','meter.s-1','meter/s' /)
  character(len=16), dimension(4), parameter :: unitnames_velkms = (/ character(len=16) :: &
  'km.s-1','km/s','kilometer.s-1','kilometer/s' /)

  ! verbosity (for debugging)
  logical, parameter :: VERBOSE = .false.

contains

! helper tools

  ! --------------------------------
  ! subroutine to check status of netcdf operation
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
  ! subroutine to check the global attributes
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
    ! remove additional remarks on units, e.g., "km (bsl)"
    if (index(val_string,"(") > 0) val_string = val_string(1:index(val_string,"(")-1)
    ! determine unit
    if (any(unitnames_m == trim(val_string))) then
      unit = 1
    else if (any(unitnames_km == trim(val_string))) then
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

  if (VERBOSE) print *

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
    print *,'Error: depth unit ',unit,' is invalid, must be either 1 == m or 2 == km'
    print *,'Please check EMC model file and attribute geospatial_vertical_units ',trim(val_string)
    stop 'EMC depth unit not recognized'
  endif

  end subroutine check_global_attributes


  ! --------------------------------
  ! subroutine to return dimension ids, names, and lengths
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
    else if (trim(dimname) == lonname_exp .or. any(lonnames == trim(dimname))) then
      lonid = dimid
      lonlen = dimlen
    else if (trim(dimname) == depname_exp .or. any(depnames == trim(dimname))) then
      depid = dimid
      deplen = dimlen
    endif
  enddo

  if (VERBOSE) print *
  if (VERBOSE) print *,'    dim latid = ', latid,' len = ',latlen
  if (VERBOSE) print *,'    dim lonid = ', lonid,' len = ',lonlen
  if (VERBOSE) print *,'    dim depid = ', depid,' len = ',deplen
  if (VERBOSE) print *

  end subroutine list_dims

  ! --------------------------------
  ! subroutine to check that the expected variable names are present
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
    else if (trim(varname) == lonname_exp .or. any(lonnames == trim(varname))) then
      varid_lon = varid
    else if (trim(varname) == depname_exp .or. any(depnames == trim(varname))) then
      varid_dep = varid
    else if (trim(varname) == vpname_exp .or. any(vpnames == trim(varname))) then
      varid_vp = varid
    else if (trim(varname) == vsname_exp .or. any(vsnames == trim(varname))) then
      varid_vs = varid
    else if (trim(varname) == rhoname_exp .or. any(rhonames == trim(varname))) then
      varid_rho = varid
    endif
  enddo

  if (VERBOSE) print *
  if (VERBOSE) print *,'    varid_lat  = ', varid_lat
  if (VERBOSE) print *,'    varid_lon  = ', varid_lon
  if (VERBOSE) print *,'    varid_dep  = ', varid_dep
  if (VERBOSE) print *
  if (VERBOSE) print *,'    varid_vp   = ', varid_vp
  if (VERBOSE) print *,'    varid_vs   = ', varid_vs
  if (VERBOSE) print *,'    varid_rho  = ', varid_rho
  if (VERBOSE) print *

  end subroutine check_varnames

  ! --------------------------------
  ! subroutine to check the dimension order of a variable
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
  latid = 0
  lonid = 0
  depid = 0
  varorderdims(:) = 0

  if (VERBOSE) print *,'  dimension ordering: '

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
      else if (trim(dimname) == lonname_exp .or. any(lonnames == trim(dimname))) then
        ! varorderdims index
        index = index + 1
        ! checks varorderdims index
        if (index > 3) stop 'Error varoderdims index exceeds 3 dimensions'
        ! Write Variable dimension order to varorderdims
        varorderdims(index) = dimids(dimid)
        islon = .true.
        lonid = index
      else if (trim(dimname) == depname_exp .or. any(depnames == trim(dimname))) then
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

  if (VERBOSE) print *
  if (VERBOSE) print *,'    Index order of dimensions for stored values: ', varorderdims
  if (VERBOSE) print *,'    lonid ',lonid,' var dim index = ',varorderdims(lonid)
  if (VERBOSE) print *,'    latid ',latid,' var dim index = ',varorderdims(latid)
  if (VERBOSE) print *,'    depid ',depid,' var dim index = ',varorderdims(depid)
  if (VERBOSE) print *

  ! checks if something found
  if (.not. islat .or. .not. islon .or. .not. isdep) then
    print *,'Error: expected dimensions lat/lon/dep are not present in variables'
    stop 'Error netcdf variable dimensions'
  endif

  end subroutine check_dimorder


  ! --------------------------------
  ! subroutine to check the attributes of a variable
  ! --------------------------------
  subroutine check_variable_attributes(ncid, varid, unit, direction, missing_val)

  implicit none
  integer, intent(in) :: ncid
  integer, intent(in) :: varid
  integer, intent(inout) :: unit,direction
  real(kind=CUSTOM_REAL), intent(inout) :: missing_val
  ! local parameters
  integer :: num_atts,i
  character(len=100) :: name,value
  real(kind=CUSTOM_REAL) :: fill_val
  logical :: has_fill_value,has_missing_value

  ! initializes
  unit        = 0
  direction   = 0
  missing_val = 9999.0_CUSTOM_REAL
  fill_val    = 9999.0_CUSTOM_REAL

  has_missing_value = .false.
  has_fill_value = .false.

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

        ! remove additional remarks on units, e.g., "km (bsl)"
        if (index(value,"(") > 0) value = value(1:index(value,"(")-1)

        ! determine unit
        if (any(unitnames_m == trim(value))) then
          unit = 1
        else if (any(unitnames_km == trim(value))) then
          unit = 2
        else if (any(unitnames_velms == trim(value))) then
          unit = 3
        else if (any(unitnames_velkms == trim(value))) then
          unit = 4
        else if (any(unitnames_gcm == trim(value))) then
          unit = 5
        else if (any(unitnames_kgcm == trim(value))) then
          unit = 6
        else if (any(unitnames_kgm == trim(value))) then
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
        has_missing_value = .true.

      else if (trim(name) == '_FillValue' .or. trim(name) == 'FillValue') then
        ! assigns value for invalid entries
        call check_status(nf90_get_att(ncid, varid, name, fill_val))
        if (VERBOSE) print *,'      fill value: ',fill_val
        has_fill_value = .true.

      else if (trim(name) == 'long_name' .or. trim(name) == '_long_name') then
        ! assigns value for invalid entries
        call check_status(nf90_get_att(ncid, varid, name, value))
        if (VERBOSE) print *,'      long name: ',trim(value)
        ! check if depth is with respect to "earth surface" or "sea level"
        if (index(value,"depth") > 0 .and. index(value,"surface") > 0) then
          ! example: long_name="depth below earth surface"
          ! depth with respect to earth surface
          EMC_depth_reference_level = 1
        else if (index(value,"depth") > 0 .and. index(value,"sea level") > 0) then
          ! example long_name="depth below sea level (bsl)"
          ! depth with respect to sea level
          EMC_depth_reference_level = 2
      endif

      endif
    enddo
  endif

  ! debug output
  if (VERBOSE) print *
  if (VERBOSE) print *,'    unit = ',unit,' direction = ',direction
  if (VERBOSE) print *,'    missing/fill value = ',missing_val,fill_val,'has missing/fill',has_missing_value,has_fill_value
  if (VERBOSE) print *

  ! sets missing value to fill value (in case only a fill value was provided)
  !
  ! note: this doesn't work for all models, as for example the AlaskaFWANT_Vs2019 model
  !       (file YangAndGao-AlaskaFWANT-Vs2019-0.0.nc) uses NaN values and not fill values
  !       to indicate missing points.
  !       we will check for NaNs, if there is no missing value provided...
  !
  if (.not. has_missing_value .and. has_fill_value) then
    missing_val = fill_val
  endif

  ! checks if something found
  if (unit == 0) then
    print *,'Error: expected attribute units is not present in variable'
    stop 'Error netcdf variable attribute'
  endif

  end subroutine check_variable_attributes

  ! --------------------------------
  ! convert NaNs to default missing_val value
  ! --------------------------------
  subroutine convert_nan_to_missing(EMC_par,missing_val)

  use, intrinsic :: ieee_arithmetic

  implicit none

  real(kind=CUSTOM_REAL), dimension(EMC_dims_nx,EMC_dims_ny,EMC_dims_nz), intent(inout) :: EMC_par
  real(kind=CUSTOM_REAL), intent(in) :: missing_val

  ! local parameters
  integer :: ix,iy,iz
  real(kind=CUSTOM_REAL) :: val

! note: isNaN() function is a GNU extension, thus not all compilers might have it.
!       using a simple val /= val check can already work.
!       here, we double-check with the IEEE intrinsic function ieee_is_nan() which should be Fortran2003 standard.

  do iz = 1,EMC_dims_nz
    do iy = 1,EMC_dims_ny
      do ix = 1,EMC_dims_nx
        ! gets value
        val = EMC_par(ix,iy,iz)

        ! checks number (isNaN)
        ! for NaN values: val /= val
        if (val /= val) EMC_par(ix,iy,iz) = missing_val

        ! double-check with ieee function
        if (ieee_support_standard(val)) then
          if (ieee_is_nan(val)) EMC_par(ix,iy,iz) = missing_val
        endif
      enddo
    enddo
  enddo

  end subroutine convert_nan_to_missing

  ! --------------------------------
  ! subroutine to scale missing density from Vp
  ! --------------------------------
  subroutine scale_Brocher_rho_from_vp()

  implicit none

  ! local parameters
  double precision :: vp,vp_p2,vp_p3,vp_p4,vp_p5,rho
  integer :: ix,iy,iz,nx,ny,nz

  ! Brocher 2005, Empirical Relations between Elastic Wavespeeds and Density in the Earth's Crust, BSSA
  ! factors from eq. (1)
  double precision,parameter :: fac1 = 1.6612d0
  double precision,parameter :: fac2 = -0.4721d0
  double precision,parameter :: fac3 = 0.0671d0
  double precision,parameter :: fac4 = -0.0043d0
  double precision,parameter :: fac5 = 0.000106d0

  ! scaling requires vp in km/s
  ! units: 1==m, 2==km, 3==m/s, 4==km/s, 5==g/cm^3, 6==kg/cm^3, 7==kg/m^3
  ! unit scaling factor
  double precision :: unit_scale

  nx = EMC_dims_nx
  ny = EMC_dims_ny
  nz = EMC_dims_nz

  ! determines scaling factor to km/s
  if (EMC_vp_unit == 3) then
    ! given in m/s -> km/s
    unit_scale = 1.d0 / 1000.d0
  else
    ! given in km/s
    unit_scale = 1.d0
  endif

  do iz = 1,nz
    do iy = 1,ny
      do ix = 1,nx
        ! Vp (in km/s)
        vp = EMC_vp(ix,iy,iz) * unit_scale
        vp_p2 = vp * vp
        vp_p3 = vp * vp_p2
        vp_p4 = vp * vp_p3
        vp_p5 = vp * vp_p4

        ! scaling relation: eq.(1)
        rho = fac1 * vp + fac2 * vp_p2 + fac3 * vp_p3 + fac4 * vp_p4 + fac5 * vp_p5

        ! Density
        EMC_rho(ix,iy,iz) = real(rho,kind=CUSTOM_REAL)
      enddo
    enddo
  enddo

  ! density scaling for rho given in g/cm^3
  EMC_rho_unit = 5

  end subroutine scale_Brocher_rho_from_vp

  ! --------------------------------
  ! subroutine to scale missing Vs from Vp
  ! --------------------------------
  subroutine scale_Brocher_vs_from_vp()

  implicit none

  ! local parameters
  double precision :: vp,vp_p2,vp_p3,vp_p4,vs
  integer :: ix,iy,iz,nx,ny,nz

  ! Brocher 2005, Empirical Relations between Elastic Wavespeeds and Density in the Earth's Crust, BSSA
  ! factors from eq. (1)
  double precision,parameter :: fac1 = 0.7858d0
  double precision,parameter :: fac2 = -1.2344d0
  double precision,parameter :: fac3 = 0.7949d0
  double precision,parameter :: fac4 = -0.1238d0
  double precision,parameter :: fac5 = 0.0064d0

  ! scaling requires vp in km/s
  ! units: 1==m, 2==km, 3==m/s, 4==km/s, 5==g/cm^3, 6==kg/cm^3, 7==kg/m^3
  ! unit scaling factor
  double precision :: unit_scale

  nx = EMC_dims_nx
  ny = EMC_dims_ny
  nz = EMC_dims_nz

  ! determines scaling factor to km/s
  if (EMC_vp_unit == 3) then
    ! given in m/s -> km/s
    unit_scale = 1.d0 / 1000.d0
  else
    ! given in km/s
    unit_scale = 1.d0
  endif

  nx = EMC_dims_nx
  ny = EMC_dims_ny
  nz = EMC_dims_nz

  do iz = 1,nz
    do iy = 1,ny
      do ix = 1,nx
        ! Vp (in km/s)
        vp = EMC_vp(ix,iy,iz) * unit_scale
        vp_p2 = vp * vp
        vp_p3 = vp * vp_p2
        vp_p4 = vp * vp_p3

        ! scaling relation: eq.(1)
        vs = fac1 + fac2 * vp + fac3 * vp_p2 + fac4 * vp_p3 + fac5 * vp_p4

        ! Density
        EMC_vs(ix,iy,iz) = real(vs,kind=CUSTOM_REAL)
      enddo
    enddo
  enddo

  ! Vs scaling for vs given in km/s
  EMC_vs_unit = 4

  ! unit scaling to convert to same unit as vp
  if (EMC_vp_unit == 3) then
    ! use same unit as vp km/s -> m/s
    EMC_vs = EMC_vs * 1000.0_CUSTOM_REAL
    EMC_vs_unit = 3
  endif

  end subroutine scale_Brocher_vs_from_vp

  ! --------------------------------
  ! subroutine to scale missing Vp from Vs
  ! --------------------------------
  subroutine scale_Brocher_vp_from_vs()

  implicit none

  ! local parameters
  double precision :: vs,vs_p2,vs_p3,vs_p4,vp
  integer :: ix,iy,iz,nx,ny,nz

  ! Brocher 2005, Empirical Relations between Elastic Wavespeeds and Density in the Earth's Crust, BSSA
  ! factors from eq. (1)
  double precision,parameter :: fac1 = 0.9409d0
  double precision,parameter :: fac2 = 2.0947d0
  double precision,parameter :: fac3 = -0.8206d0
  double precision,parameter :: fac4 = 0.2683d0
  double precision,parameter :: fac5 = -0.0251d0

  ! scaling requires vs in km/s
  ! units: 1==m, 2==km, 3==m/s, 4==km/s, 5==g/cm^3, 6==kg/cm^3, 7==kg/m^3
  ! unit scaling factor
  double precision :: unit_scale

  nx = EMC_dims_nx
  ny = EMC_dims_ny
  nz = EMC_dims_nz

  ! determines scaling factor to km/s
  if (EMC_vs_unit == 3) then
    ! given in m/s -> km/s
    unit_scale = 1.d0 / 1000.d0
  else
    ! given in km/s
    unit_scale = 1.d0
  endif

  nx = EMC_dims_nx
  ny = EMC_dims_ny
  nz = EMC_dims_nz

  do iz = 1,nz
    do iy = 1,ny
      do ix = 1,nx
        ! Vs (in km/s)
        vs = EMC_vs(ix,iy,iz) * unit_scale
        vs_p2 = vs * vs
        vs_p3 = vs * vs_p2
        vs_p4 = vs * vs_p3

        ! scaling relation: eq.(1)
        vp = fac1 + fac2 * vs + fac3 * vs_p2 + fac4 * vs_p3 + fac5 * vs_p4

        ! Vp
        EMC_vp(ix,iy,iz) = real(vp,kind=CUSTOM_REAL)
      enddo
    enddo
  enddo

  ! Vp scaling for vp given in km/s
  EMC_vp_unit = 4

  ! unit scaling to convert to same unit as vs
  if (EMC_vs_unit == 3) then
    ! use same unit as vp km/s -> m/s
    EMC_vp = EMC_vp * 1000.0_CUSTOM_REAL
    EMC_vp_unit = 3
  endif

  end subroutine scale_Brocher_vp_from_vs

end module model_emc_par

!-------------------------------------------------------------------------------------------
!
! helper routines
!
!-------------------------------------------------------------------------------------------

  subroutine reorder_EMC_model_arrays(lonid,latid,depid,varorderdims)

  use model_emc_par

  implicit none

  integer, intent(in) :: latid, lonid, depid       ! dimension ids
  integer, intent(in) :: varorderdims(3)           ! variable dimension order

  ! local parameters
  integer :: ilon,ilat,idep,ier
  integer :: ix,iy,iz
  integer :: dimindex(3)
  ! temporary arrays for re-ordering
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: tmp_array
  logical, dimension(:,:,:), allocatable :: tmp_mask

  ! note: arrays EMC_vp(:,:,:),.. might use different index orders for different EMC model files.
  !       For example, Alaska file uses nx==lonlen, ny==latlen, nz==deplen,
  !       but other EMC model files might have a different order.
  !       This gets determines by dimlens() and varorderdims().
  !
  ! user output
  write(IMAIN,*) '  reordering:'
  write(IMAIN,*) '    current indexing: lon/lat/depth ids = ',lonid,latid,depid
  write(IMAIN,*) '                                var ids = ',varorderdims(lonid),varorderdims(latid),varorderdims(depid)
  write(IMAIN,*) '                dimensions nx / ny / nz = ',EMC_dims_nx,EMC_dims_ny,EMC_dims_nz
  call flush_IMAIN()

  ! checks if longitude needs re-ordering as well
  if (EMC_dlon < 0.0) then
    write(IMAIN,*) '                lon ordering  : decreases ',EMC_dlon
    call flush_IMAIN()
  endif

  ! checks if latitude needs re-ordering as well
  if (EMC_dlat < 0.0) then
    write(IMAIN,*) '                lat ordering  : decreases ',EMC_dlat
    call flush_IMAIN()
  endif

  ! checks if depth needs re-ordering as well
  if (EMC_ddep < 0.0) then
    write(IMAIN,*) '                depth ordering: decreases ',EMC_ddep
    call flush_IMAIN()
  endif

  ! allocates temporary array
  allocate(tmp_array(EMC_lonlen,EMC_latlen,EMC_deplen),stat=ier)
  if (ier /= 0) stop 'Error allocating temporary array'
  tmp_array(:,:,:) = 0.0_CUSTOM_REAL

  ! Vp array
  do idep = 1,EMC_deplen
    do ilat = 1,EMC_latlen
      do ilon = 1,EMC_lonlen
        ! gets array indexing
        if (EMC_dlon < 0.0) then
          dimindex(varorderdims(lonid)) = EMC_lonlen - ilon + 1
        else
          dimindex(varorderdims(lonid)) = ilon
        endif
        if (EMC_dlat < 0.0) then
          dimindex(varorderdims(latid)) = EMC_latlen - ilat + 1
        else
          dimindex(varorderdims(latid)) = ilat
        endif
        if (EMC_ddep < 0.0) then
          dimindex(varorderdims(depid)) = EMC_deplen - idep + 1
        else
          dimindex(varorderdims(depid)) = idep
        endif
        ix = dimindex(varorderdims(1))
        iy = dimindex(varorderdims(2))
        iz = dimindex(varorderdims(3))
        ! saves in re-order temporary array
        tmp_array(ilon,ilat,idep) = EMC_vp(ix,iy,iz)
      enddo
    enddo
  enddo
  ! re-allocate with new order
  deallocate(EMC_vp)
  allocate(EMC_vp(EMC_lonlen,EMC_latlen,EMC_deplen),stat=ier)
  if (ier /= 0) stop 'Error allocating reordered vp array'
  EMC_vp(:,:,:) = tmp_array(:,:,:)

  ! Vs array
  tmp_array(:,:,:) = 0.0_CUSTOM_REAL
  do idep = 1,EMC_deplen
    do ilat = 1,EMC_latlen
      do ilon = 1,EMC_lonlen
        ! gets array indexing
        if (EMC_dlon < 0.0) then
          dimindex(varorderdims(lonid)) = EMC_lonlen - ilon + 1
        else
          dimindex(varorderdims(lonid)) = ilon
        endif
        if (EMC_dlat < 0.0) then
          dimindex(varorderdims(latid)) = EMC_latlen - ilat + 1
        else
          dimindex(varorderdims(latid)) = ilat
        endif
        if (EMC_ddep < 0.0) then
          dimindex(varorderdims(depid)) = EMC_deplen - idep + 1
        else
          dimindex(varorderdims(depid)) = idep
        endif
        ix = dimindex(varorderdims(1))
        iy = dimindex(varorderdims(2))
        iz = dimindex(varorderdims(3))
        ! saves in re-order temporary array
        tmp_array(ilon,ilat,idep) = EMC_vs(ix,iy,iz)
      enddo
    enddo
  enddo
  ! re-allocate with new order
  deallocate(EMC_vs)
  allocate(EMC_vs(EMC_lonlen,EMC_latlen,EMC_deplen),stat=ier)
  if (ier /= 0) stop 'Error allocating reordered vs array'
  EMC_vs(:,:,:) = tmp_array(:,:,:)

  ! density array
  tmp_array(:,:,:) = 0.0_CUSTOM_REAL
  do idep = 1,EMC_deplen
    do ilat = 1,EMC_latlen
      do ilon = 1,EMC_lonlen
        ! gets array indexing
        if (EMC_dlon < 0.0) then
          dimindex(varorderdims(lonid)) = EMC_lonlen - ilon + 1
        else
          dimindex(varorderdims(lonid)) = ilon
        endif
        if (EMC_dlat < 0.0) then
          dimindex(varorderdims(latid)) = EMC_latlen - ilat + 1
        else
          dimindex(varorderdims(latid)) = ilat
        endif
        if (EMC_ddep < 0.0) then
          dimindex(varorderdims(depid)) = EMC_deplen - idep + 1
        else
          dimindex(varorderdims(depid)) = idep
        endif
        ix = dimindex(varorderdims(1))
        iy = dimindex(varorderdims(2))
        iz = dimindex(varorderdims(3))
        ! saves in re-order temporary array
        tmp_array(ilon,ilat,idep) = EMC_rho(ix,iy,iz)
      enddo
    enddo
  enddo
  ! re-allocate with new order
  deallocate(EMC_rho)
  allocate(EMC_rho(EMC_lonlen,EMC_latlen,EMC_deplen),stat=ier)
  if (ier /= 0) stop 'Error allocating reordered vs array'
  EMC_rho(:,:,:) = tmp_array(:,:,:)

  ! re-order mask
  ! (mask is .true. for missing value)
  allocate(tmp_mask(EMC_lonlen,EMC_latlen,EMC_deplen),stat=ier)
  if (ier /= 0) stop 'Error allocating temporary mask array'
  tmp_mask(:,:,:) = .false.
  do idep = 1,EMC_deplen
    do ilat = 1,EMC_latlen
      do ilon = 1,EMC_lonlen
        ! gets array indexing
        if (EMC_dlon < 0.0) then
          dimindex(varorderdims(lonid)) = EMC_lonlen - ilon + 1
        else
          dimindex(varorderdims(lonid)) = ilon
        endif
        if (EMC_dlat < 0.0) then
          dimindex(varorderdims(latid)) = EMC_latlen - ilat + 1
        else
          dimindex(varorderdims(latid)) = ilat
        endif
        if (EMC_ddep < 0.0) then
          dimindex(varorderdims(depid)) = EMC_deplen - idep + 1
        else
          dimindex(varorderdims(depid)) = idep
        endif
        ix = dimindex(varorderdims(1))
        iy = dimindex(varorderdims(2))
        iz = dimindex(varorderdims(3))
        ! saves in re-order temporary array
        tmp_mask(ilon,ilat,idep) = EMC_mask(ix,iy,iz)
      enddo
    enddo
  enddo
  ! re-allocate with new order
  deallocate(EMC_mask)
  allocate(EMC_mask(EMC_lonlen,EMC_latlen,EMC_deplen),stat=ier)
  if (ier /= 0) stop 'Error allocating reordered mask array'
  EMC_mask(:,:,:) = tmp_mask(:,:,:)

  ! re-orders lon array
  if (EMC_dlon < 0.0) then
    tmp_array(:,:,:) = 0.0_CUSTOM_REAL
    do ilon = 1,EMC_lonlen
      tmp_array(ilon,1,1) = EMC_lon(EMC_lonlen - ilon + 1)
    enddo
    EMC_lon(:) = tmp_array(:,1,1)
    EMC_dlon = - EMC_dlon
  endif

  ! re-orders lat array
  if (EMC_dlat < 0.0) then
    tmp_array(:,:,:) = 0.0_CUSTOM_REAL
    do ilat = 1,EMC_latlen
      tmp_array(1,ilat,1) = EMC_lat(EMC_latlen - ilat + 1)
    enddo
    EMC_lat(:) = tmp_array(1,:,1)
    EMC_dlat = - EMC_dlat
  endif

  ! re-orders depth array
  if (EMC_ddep < 0.0) then
    tmp_array(:,:,:) = 0.0_CUSTOM_REAL
    do idep = 1,EMC_deplen
      tmp_array(1,1,idep) = EMC_dep(EMC_deplen - idep + 1)
    enddo
    EMC_dep(:) = tmp_array(1,1,:)
    EMC_ddep = - EMC_ddep
  endif

  ! free temporary array
  deallocate(tmp_array)
  deallocate(tmp_mask)

  ! re-set dimensions
  EMC_dims_nx = EMC_lonlen
  EMC_dims_ny = EMC_latlen
  EMC_dims_nz = EMC_deplen

  ! user output
  write(IMAIN,*)
  write(IMAIN,*) '    new     indexing: lon/lat/depth ids = ',1,2,3
  write(IMAIN,*) '                dimensions nx / ny / nz = ',EMC_dims_nx,EMC_dims_ny,EMC_dims_nz
  write(IMAIN,*) '                lon ordering  :    dlon = ',EMC_dlon
  write(IMAIN,*) '                lat ordering  :    dlat = ',EMC_dlat
  write(IMAIN,*) '                depth ordering:    ddep = ',EMC_ddep
  write(IMAIN,*)
  call flush_IMAIN()

  end subroutine reorder_EMC_model_arrays

!
!-------------------------------------------------------------------------------------------
!

  subroutine fill_EMC_surface_values()

  use model_emc_par

  implicit none

  ! local parameters
  integer :: ilon,ilat,idep,icount
  integer :: surface_index
  real(kind=CUSTOM_REAL) :: val_vp,val_vs,val_rho

  ! to avoid taking model values defined for air or oceans, we search for the top-most index that defines the
  ! velocities for the solid surface and extend those solid values upwards

  ! check if something to do
  if (count(abs(EMC_vs(:,:,:)) < TINYVAL) == 0) return

  ! user output
  write(IMAIN,*) '  filling fluid surface values w/ solid:'
  call flush_IMAIN()

  ! counter
  icount = 0

  ! for each lat/lon grid position
  do ilat = 1,EMC_latlen
    do ilon = 1,EMC_lonlen
      ! gets first index where Vs > 0
      ! (assumes that depth array starts at top and moves down to earth center)
      surface_index = 0
      do idep = 1,EMC_deplen
        ! checks if point has valid entry
        if (.not. EMC_mask(ilon,ilat,idep) .and. surface_index == 0) then
          ! checks if point is for solids (vs > 0)
          if (EMC_vs(ilon,ilat,idep) > TINYVAL) then
            ! store index
            surface_index = idep
            ! exit depth loop
            exit
          endif
        endif
      enddo

      ! checks missed array values (in case some lat/lon are all missing values)
      if (surface_index == 0) surface_index = 1

      ! extends surface velocities upwards
      ! EMC models might have values for air/water on top - for such points above solid surface elevation,
      ! we extend the values corresponding to the solid surface (no oceans, air, etc.)
      ! solid values
      val_vp = EMC_vp(ilon,ilat,surface_index)
      val_vs = EMC_vs(ilon,ilat,surface_index)
      val_rho = EMC_rho(ilon,ilat,surface_index)
      ! extend values upwards
      do idep = 1,surface_index
        EMC_vp(ilon,ilat,idep) = val_vp
        EMC_vs(ilon,ilat,idep) = val_vs
        EMC_rho(ilon,ilat,idep) = val_rho
      enddo

      ! updates entries per array
      icount = icount + (surface_index - 1)
    enddo
  enddo

  ! user output
  write(IMAIN,*) '    updated surface entries = ',icount
  if (icount == 0) write(IMAIN,*) '    no fluid (vs == 0) entries found'
  write(IMAIN,*)
  call flush_IMAIN()

  end subroutine fill_EMC_surface_values


!
!-------------------------------------------------------------------------------------------
!

  subroutine determine_EMC_depth_average()

  use constants, only: IOUT
  use shared_parameters, only: OUTPUT_FILES

  use model_emc_par

  implicit none
  ! local parameters
  integer :: idep,ier
  character (len=MAX_STRING_LEN) :: filename
  ! spline interpolation
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: spline_b_vp,spline_c_vp,spline_d_vp
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: spline_b_vs,spline_c_vs,spline_d_vs
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: spline_b_rho,spline_c_rho,spline_d_rho
  real(kind=CUSTOM_REAL) :: val_interp_vp,val_interp_vs,val_interp_rho
  real(kind=CUSTOM_REAL) :: z_depth
  integer :: i
  ! for spline interpolation of depth profiles
  logical, parameter :: OUTPUT_SPLINE_INTERPOLATED_DEPTH = .false.

  ! (based on defined grid point values only, ignoring missing point values)
  ! user output
  write(IMAIN,*) '  EMC model 1D depth average:'
  call flush_IMAIN()

  ! allocates 1D average arrays
  allocate(EMC_avg_rho(EMC_deplen),EMC_avg_vp(EMC_deplen),EMC_avg_vs(EMC_deplen),stat=ier)
  if (ier /= 0) stop 'Error allocating 1d average arrays'
  EMC_avg_rho(:) = 0.0; EMC_avg_vp(:) = 0.0; EMC_avg_vs(:) = 0.0

  ! Vp
  call get_average_model(EMC_vp,EMC_avg_vp)
  ! Vs
  call get_average_model(EMC_vs,EMC_avg_vs)
  ! Density
  call get_average_model(EMC_rho,EMC_avg_rho)

  ! save output as file
  filename = trim(OUTPUT_FILES) // '/EMC_model_1D_depth_average.dat'
  open(unit=IOUT,file=trim(filename),status='unknown',iostat=ier)
  if (ier /= 0) call exit_mpi(0,'Error opening EMC_model_1D_depth_average.dat file')
  write(IOUT,*) '# EMC model - 1D depth average'
  write(IOUT,*) '#'
  write(IOUT,*) '# model id    : ',trim(EMC_model_id)
  write(IOUT,*) '#       title : ',trim(EMC_model_title)
  write(IOUT,*) '#       range : latitude  min/max = ',EMC_lat_min,'/',EMC_lat_max
  write(IOUT,*) '#               longitude min/max = ',EMC_lon_min,'/',EMC_lon_max
  write(IOUT,*) '#               depth     min/max = ',EMC_dep_min,'/',EMC_dep_max,'(km)'
  write(IOUT,*) '#'
  write(IOUT,*) '# 1D depth average'
  write(IOUT,*) '# format:'
  write(IOUT,*) '#depth (km) #rho (kg/m^3)  #vp (m/s)  #vs (m/s)'
  do idep = 1,EMC_deplen
    write(IOUT,*) EMC_dep(idep),EMC_avg_rho(idep),EMC_avg_vp(idep),EMC_avg_vs(idep)
  enddo
  close(IOUT)

  ! user output
  write(IMAIN,*) '    written to: ',trim(filename)
  write(IMAIN,*)
  call flush_IMAIN()

  ! spline evaluations
  if (OUTPUT_SPLINE_INTERPOLATED_DEPTH) then
    ! allocates 1D average arrays
    allocate(spline_b_vp(EMC_deplen),spline_c_vp(EMC_deplen),spline_d_vp(EMC_deplen),stat=ier)
    if (ier /= 0) stop 'Error allocating spline vp arrays'
    spline_b_vp(:) = 0.0; spline_c_vp(:) = 0.0; spline_d_vp(:) = 0.0

    allocate(spline_b_vs(EMC_deplen),spline_c_vs(EMC_deplen),spline_d_vs(EMC_deplen),stat=ier)
    if (ier /= 0) stop 'Error allocating spline vs arrays'
    spline_b_vs(:) = 0.0; spline_c_vs(:) = 0.0; spline_d_vs(:) = 0.0

    allocate(spline_b_rho(EMC_deplen),spline_c_rho(EMC_deplen),spline_d_rho(EMC_deplen),stat=ier)
    if (ier /= 0) stop 'Error allocating spline rho arrays'
    spline_b_rho(:) = 0.0; spline_c_rho(:) = 0.0; spline_d_rho(:) = 0.0

    call EMC_setup_spline_coeffs(EMC_deplen,EMC_dep,EMC_avg_vp,spline_b_vp,spline_c_vp,spline_d_vp)
    call EMC_setup_spline_coeffs(EMC_deplen,EMC_dep,EMC_avg_vs,spline_b_vs,spline_c_vs,spline_d_vs)
    call EMC_setup_spline_coeffs(EMC_deplen,EMC_dep,EMC_avg_rho,spline_b_rho,spline_c_rho,spline_d_rho)

    ! save output as file
    filename = trim(OUTPUT_FILES) // '/EMC_model_1D_depth_average.spline-interpolated.dat'
    open(unit=IOUT,file=trim(filename),status='unknown',iostat=ier)
    write(IOUT,*) '# 1D depth average - interpolated'
    write(IOUT,*) '#depth (km) #rho (kg/m^3)  #vp (m/s)  #vs (m/s)'

    ! interpolate
    do i = 0,500
      z_depth = i * (EMC_dep_max / 500.0)

      ! Perform cubic B-spline interpolation
      call EMC_eval_spline(EMC_deplen,EMC_dep,EMC_avg_vp,spline_b_vp,spline_c_vp,spline_d_vp,z_depth,val_interp_vp)
      call EMC_eval_spline(EMC_deplen,EMC_dep,EMC_avg_vs,spline_b_vs,spline_c_vs,spline_d_vs,z_depth,val_interp_vs)
      call EMC_eval_spline(EMC_deplen,EMC_dep,EMC_avg_rho,spline_b_rho,spline_c_rho,spline_d_rho,z_depth,val_interp_rho)

      write(IOUT,*) z_depth,val_interp_rho,val_interp_vp,val_interp_vs
    enddo
    close(IOUT)
  endif

contains

  subroutine get_average_model(EMC_par,EMC_avg)

  use model_emc_par, only: EMC_dims_nx,EMC_dims_ny,EMC_dims_nz,EMC_mask

  implicit none
  real(kind=CUSTOM_REAL), dimension(EMC_dims_nx,EMC_dims_ny,EMC_dims_nz), intent(in) :: EMC_par
  real(kind=CUSTOM_REAL), dimension(EMC_deplen), intent(inout) :: EMC_avg
  ! local variables
  integer :: idep,ilat,ilon,icount
  real(kind=CUSTOM_REAL) :: val_avg,val

  ! initializes average value
  EMC_avg(:) = 0.0_CUSTOM_REAL

  ! loops over each depth and determines the average value for valid points
  do idep = 1,EMC_deplen
    ! initializes depth average
    icount = 0
    val_avg = 0.0_CUSTOM_REAL

    ! loops over horizontal slice
    do ilat = 1,EMC_latlen
      do ilon = 1,EMC_lonlen
        ! checks if point has valid entry
        ! (mask is .true. for missing points)
        if (.not. EMC_mask(ilon,ilat,idep)) then
          ! gets model parameter value (vp, vs or rho)
          val = EMC_par(ilon,ilat,idep)
          ! checks if point is for solids (vs > 0)
          if (val > TINYVAL) then
            ! adds value
            icount = icount + 1
            val_avg = val_avg + val
          endif
        endif
      enddo
    enddo

    ! stores average value
    if (icount > 0) EMC_avg(idep) = val_avg / icount

    !debug
    !print *,'debug: average ',idep,EMC_avg(idep),' - counts',icount
  enddo

  end subroutine get_average_model

  end subroutine determine_EMC_depth_average

!
!-------------------------------------------------------------------------------------------
!

  subroutine EMC_setup_spline_coeffs(n,x,y,spline_b,spline_c,spline_d)

  use constants, only: CUSTOM_REAL

  implicit none
  integer, intent(in) :: n
  real(kind=CUSTOM_REAL), dimension(n), intent(in) :: x, y  ! x position, y values
  real(kind=CUSTOM_REAL), dimension(n), intent(inout) :: spline_b,spline_c,spline_d

  ! local parameters
  integer :: i, j
  real(kind=CUSTOM_REAL) :: a(n), b(n), c(n)
  real(kind=CUSTOM_REAL) :: h(n-1), d(n), l(n), mu(n), z(n)

  ! initializes
  a(:) = 0.0
  b(:) = 0.0
  c(:) = 0.0
  d(:) = 0.0
  l(:) = 0.0
  mu(:) = 0.0
  z(:) = 0.0

  ! Compute h[i] = x[i+1] - x[i] and initialize arrays
  do i = 1, n-1
    h(i) = x(i+1) - x(i)
  enddo

  ! Compute a[i] = (3/h[i]) * (a[i+1] - a[i]) - (3/h[i-1]) * (a[i] - a[i-1])
  do i = 2, n-1
    a(i) = (3.0 / h(i)) * (y(i+1) - y(i)) - (3.0 / h(i-1)) * (y(i) - y(i-1))
  enddo

  ! Solve for c[i] using Thomas algorithm for tridiagonal systems
  l(1) = 2.0 * h(1)
  mu(1) = 0.5
  z(1) = a(1) / l(1)
  do i = 2, n-1
    l(i) = 2.0 * (x(i+1) - x(i-1)) - h(i-1) * mu(i-1)
    mu(i) = h(i) / l(i)
    z(i) = (a(i) - h(i-1) * z(i-1)) / l(i)
  enddo
  l(n) = h(n-1) * (2.0 - mu(n-1))
  z(n) = (a(n) - h(n-1) * z(n-1)) / l(n)
  c(n) = z(n)
  do j = n-1, 1, -1
    c(j) = z(j) - mu(j) * c(j+1)
    b(j) = (y(j+1) - y(j)) / h(j) - h(j) * (c(j+1) + 2.0 * c(j)) / 3.0
    d(j) = (c(j+1) - c(j)) / (3.0 * h(j))
  enddo

  ! return coefficients
  spline_b(:) = b(:)
  spline_c(:) = c(:)
  spline_d(:) = d(:)

  end subroutine EMC_setup_spline_coeffs

!
!-------------------------------------------------------------------------------------------
!

  subroutine EMC_eval_spline(n,x,y,spline_b,spline_c,spline_d,x_target,y_interp)

  use constants, only: CUSTOM_REAL

  implicit none
  integer, intent(in) :: n
  real(kind=CUSTOM_REAL), dimension(n), intent(in) :: x,y
  real(kind=CUSTOM_REAL), dimension(n), intent(in) :: spline_b,spline_c,spline_d
  real(kind=CUSTOM_REAL), intent(in) :: x_target
  real(kind=CUSTOM_REAL), intent(inout) :: y_interp
  ! local parameters
  integer :: i

  ! Find the segment that x_target is in and evaluate the spline there
  i = 1
  do while (i < n)
    if (x_target >= x(i) .and. x_target <= x(i+1)) exit
    i = i + 1
  enddo

  ! spline interpolation
  y_interp = y(i) + spline_b(i) * (x_target - x(i)) + spline_c(i) * (x_target - x(i))**2 + spline_d(i) * (x_target - x(i))**3

  end subroutine EMC_eval_spline

!
!-------------------------------------------------------------------------------------------
!

  subroutine fill_EMC_missing_values_average()

! replaces missing values with 1D depth average values

  use model_emc_par

  implicit none
  ! local parameters
  integer :: ilat,ilon,idep,icount
  real(kind=CUSTOM_REAL) :: val_vs,val_vp,val_rho

  ! checks if anything to do
  if (count(EMC_mask(:,:,:) .eqv. .true.) == 0) return

  ! user output
  write(IMAIN,*) '  filling:'
  write(IMAIN,*) '    using 1D depth average for missing values'
  call flush_IMAIN()

  ! counter
  icount = 0

  ! loops over each depth and replaces missing values with the average 1D value
  do idep = 1,EMC_deplen
    ! gets depth average values
    val_vp = EMC_avg_vp(idep)
    val_vs = EMC_avg_vs(idep)
    val_rho = EMC_avg_rho(idep)

    ! loops over horizontal slice
    do ilat = 1,EMC_latlen
      do ilon = 1,EMC_lonlen
        ! replace missing value
        if (EMC_mask(ilon,ilat,idep)) then
          EMC_vp(ilon,ilat,idep) = val_vp
          EMC_vs(ilon,ilat,idep) = val_vs
          EMC_rho(ilon,ilat,idep) = val_rho
          ! re-set mask
          EMC_mask(ilon,ilat,idep) = .false.
          ! counter
          icount = icount + 1
        endif
      enddo
    enddo
  enddo

  ! checks that all values have been set
  if (count(EMC_mask(:,:,:) .eqv. .true.) > 0) stop 'Error filling missing values w/ average'

  ! user output
  write(IMAIN,*) '    number of updated values = ',icount
  write(IMAIN,*) '    vp  min/max = ',minval(EMC_vp),'/',maxval(EMC_vp),'(m/s)'
  write(IMAIN,*) '    vs  min/max = ',minval(EMC_vs),'/',maxval(EMC_vs),'(m/s)'
  write(IMAIN,*) '    rho min/max = ',minval(EMC_rho),'/',maxval(EMC_rho),'(kg/m^3)'
  write(IMAIN,*)
  call flush_IMAIN()

  end subroutine fill_EMC_missing_values_average

!
!-------------------------------------------------------------------------------------------
!

  subroutine fill_EMC_missing_values_interpolated()

! replaces missing values with interpolated values from neighboring grid points

  use model_emc_par

  implicit none

  ! local parameters
  integer :: ilat,ilon,idep,iupdated,ier
  integer :: Nx,Ny
  real(kind=CUSTOM_REAL) :: val_avg_vp,val_avg_vs,val_avg_rho
  double precision :: vp_interp,vs_interp,rho_interp
  logical, dimension(:,:,:), allocatable :: tmp_mask

  ! grid search number of grid steps (Shepard & nearest neighbor)
  integer :: NUMSTEPS_IX,NUMSTEPS_IY
  integer, parameter :: NUMSTEPS = 8

  ! tapering
  double precision :: taper_dist_maximum

  ! checks if anything to do
  if (count(EMC_mask(:,:,:) .eqv. .true.) == 0) return

  ! user output
  write(IMAIN,*) '  filling:'
  write(IMAIN,*) '    using interpolated closest model values for missing values'
  write(IMAIN,*) '    interpolation method : ',INTERPOLATION_METHOD,'(1 == Shepard/2 == nearest/3 == bilinear)'
  write(IMAIN,*)
  call flush_IMAIN()

  ! allocate temporary mask copy
  allocate(tmp_mask(EMC_dims_nx,EMC_dims_ny,EMC_dims_nz),stat=ier)
  if (ier /= 0) stop 'Error allocating tmp_mask array'
  ! copy mask
  tmp_mask(:,:,:) = EMC_mask(:,:,:)

  ! counter
  iupdated = 0

  Nx = EMC_lonlen
  Ny = EMC_latlen

  ! determines number of grid steps for Shepard & nearest neighbor and different dlon,dlat sizes
  NUMSTEPS_IX = nint( NUMSTEPS * max(EMC_dlon,EMC_dlat) / EMC_dlon )
  NUMSTEPS_IY = nint( NUMSTEPS * max(EMC_dlon,EMC_dlat) / EMC_dlat )

  ! maximum distance for tapering
  taper_dist_maximum = NUMSTEPS * max(EMC_dlon,EMC_dlat)

  ! loops over each depth and replaces missing values with values interpolated from neighboring grid points
  do idep = 1,EMC_deplen
    ! gets depth average values
    val_avg_vp = EMC_avg_vp(idep)
    val_avg_vs = EMC_avg_vs(idep)
    val_avg_rho = EMC_avg_rho(idep)

    ! loops over horizontal slice
    do ilat = 1,EMC_latlen
      do ilon = 1,EMC_lonlen
        ! replace missing value
        if (EMC_mask(ilon,ilat,idep)) then
          ! bilinear interpolation
          select case (INTERPOLATION_METHOD)
          case (1)
            ! inverse distance weighted (IDL) interpolation (Shepard's method)
            call do_shepard_interpolation(vp_interp,vs_interp,rho_interp)
          case (2)
            ! nearest neighbor
            call do_nearest_interpolation(vp_interp,vs_interp,rho_interp)
          case (3)
            ! bilinear
            call do_bilinear_interpolation(vp_interp,vs_interp,rho_interp)
          case default
            stop 'Interpolation method not implemented for filling missing values'
          end select

          ! stores interpolated value
          EMC_vp(ilon,ilat,idep) = real(vp_interp,kind=CUSTOM_REAL)
          EMC_vs(ilon,ilat,idep) = real(vs_interp,kind=CUSTOM_REAL)
          EMC_rho(ilon,ilat,idep) = real(rho_interp,kind=CUSTOM_REAL)

          ! update mask flag
          tmp_mask(ilon,ilat,idep) = .false.

          ! counter
          iupdated = iupdated + 1
        endif
      enddo
    enddo
  enddo

  ! update main mask, should be all .false. now
  EMC_mask(:,:,:) = tmp_mask(:,:,:)

  ! checks that all values have been set
  if (count(EMC_mask(:,:,:) .eqv. .true.) > 0) stop 'Error filling missing values w/ average'

  ! user output
  write(IMAIN,*) '    number of updated values = ',iupdated
  write(IMAIN,*) '    vp  min/max = ',minval(EMC_vp),'/',maxval(EMC_vp),'(m/s)'
  write(IMAIN,*) '    vs  min/max = ',minval(EMC_vs),'/',maxval(EMC_vs),'(m/s)'
  write(IMAIN,*) '    rho min/max = ',minval(EMC_rho),'/',maxval(EMC_rho),'(kg/m^3)'
  write(IMAIN,*)
  call flush_IMAIN()

  ! free temporary array
  deallocate(tmp_mask)

contains

  function get_next_index_lon(ilon,ilat,idep,idirection) result (ix)

  implicit none
  integer, intent(in) :: ilon,ilat,idep,idirection
  ! return result
  integer :: ix
  ! local parameters

  ! find defined point
  ! (mask is .true. for missing value)
  ix = ilon

  ! direction -1 == search to the left/smaller, +1 == search to the right/bigger
  select case(idirection)
  case (-1)
    ! looks for smaller index that is defined
    do while (EMC_mask(ix,ilat,idep))
      ix = ix - 1
      if (ix < 1) exit
    enddo

  case (1)
    ! looks for bigger index that is defined
    do while (EMC_mask(ix,ilat,idep))
      ix = ix + 1
      if (ix > EMC_lonlen) exit
    enddo
  end select

  end function get_next_index_lon

  !------------------------------------

  function get_next_index_lat(ilon,ilat,idep,idirection) result (iy)

  implicit none
  integer, intent(in) :: ilon,ilat,idep,idirection
  ! return result
  integer :: iy
  ! local parameters

  ! find defined point
  ! (mask is .true. for missing value)
  iy = ilat

  ! direction -1 == search to the left/smaller, +1 == search to the right/bigger
  select case(idirection)
  case (-1)
    ! looks for smaller index that is defined
    do while (EMC_mask(ilon,iy,idep))
      iy = iy - 1
      if (iy < 1) exit
    enddo

  case (1)
    ! looks for bigger index that is defined
    do while (EMC_mask(ilon,iy,idep))
      iy = iy + 1
      if (iy > EMC_latlen) exit
    enddo
  end select

  end function get_next_index_lat

  !------------------------------------

  subroutine do_nearest_interpolation(vp_interp,vs_interp,rho_interp)

  use constants, only: HUGEVAL

  implicit none

  double precision, intent(inout) :: vp_interp,vs_interp,rho_interp

  ! local parameters
  double precision :: lon_ref,lat_ref,lon,lat,dlon,dlat
  double precision :: dist_sq
  integer :: i,j,ix,iy
  ! nearest neighbor
  integer :: ilon_n,ilat_n
  double precision :: dist_n
  ! tapering
  double precision :: taper_val
  double precision, external :: cosine_taper

  ! reference grid point position
  lon_ref = EMC_lon(ilon)
  lat_ref = EMC_lat(ilat)

  dist_n = +HUGEVAL
  ilon_n = 0
  ilat_n = 0

  ! finds nearest point within a grid search
  do j = -NUMSTEPS_IY,NUMSTEPS_IY
    do i = -NUMSTEPS_IX,NUMSTEPS_IX
      ! grid index
      ix = ilon + i
      iy = ilat + j

      ! checks bounds
      if (ix < 1) cycle
      if (ix > Nx) cycle
      if (iy < 1) cycle
      if (iy > Ny) cycle

      ! skip missing points
      if (EMC_mask(ix,iy,idep)) cycle

      ! defined point location
      lon = EMC_lon(ix)
      lat = EMC_lat(iy)

      ! distance - as simple as it gets in degrees
      ! no Haversine function (that could be used for more accurate results, but more costly)
      dlon = lon_ref - lon
      dlat = lat_ref - lat
      dist_sq = dlon*dlon + dlat*dlat
      if (dist_sq < dist_n) then
        dist_n = dist_sq
        ilon_n = ix
        ilat_n = iy
      endif
    enddo
  enddo

  !debug
  !print *,'debug: nearest interpolation: dist_n ',dist_n,' ilon_n/ilat_n ',ilon_n,ilat_n,'Nx/Ny',Nx,Ny
  !print *,'       reference lon/lat',lon_ref,lat_ref
  !print *,'       EMC lon min/max',EMC_lon_min,EMC_lon_max
  !print *,'       EMC lat min/max',EMC_lat_min,EMC_lat_max

  ! checks if we found a point
  if (ilon_n == 0 .or. ilat_n == 0) then
    ! no near point found
    ! return background/average value
    vp_interp = val_avg_vp
    vs_interp = val_avg_vs
    rho_interp = val_avg_rho
  else
    ! computes taper
    dist_n = sqrt(dist_n)
    taper_val = cosine_taper(dist_n, taper_dist_maximum)

    ! nearest point value w/ tapering
    vp_interp = EMC_vp(ilon_n,ilat_n,idep) * taper_val + (1.d0 - taper_val) * val_avg_vp
    vs_interp = EMC_vs(ilon_n,ilat_n,idep) * taper_val + (1.d0 - taper_val) * val_avg_vs
    rho_interp = EMC_rho(ilon_n,ilat_n,idep) * taper_val + (1.d0 - taper_val) * val_avg_rho
  endif

  end subroutine do_nearest_interpolation

  !------------------------------------

  subroutine do_bilinear_interpolation(vp_interp,vs_interp,rho_interp)

  implicit none

  double precision, intent(inout) :: vp_interp,vs_interp,rho_interp

  ! local parameters
  double precision :: lon,lat
  double precision :: dx,dy
  double precision :: gamma_interp_x,gamma_interp_y
  double precision :: val1_par(3),val2_par(3),val3_par(3),val4_par(3)
  integer :: ix,iy,ixplus,iyplus

  ! to use average/background velocities for "virtual" points on boarders
  logical, parameter :: FILL_MISSING_BORDER_POINTS_WITH_AVERAGE_VALUES = .true.

  ! get lon index from next grid point that is defined
  ! left
  ix = get_next_index_lon(ilon,ilat,idep,-1)
  ! right
  ixplus = get_next_index_lon(ilon,ilat,idep,+1)

  ! get lat index from next grid point that is defined
  ! down
  iy = get_next_index_lat(ilon,ilat,idep,-1)
  ! up
  iyplus = get_next_index_lat(ilon,ilat,idep,+1)

  ! sets interpolation values
  ! longitude
  if (ix < 1 .and. ixplus > Nx) then
    ! no valid point on this lon
    ! puts "virtual" points to border
    ix = 1
    ixplus = Nx
    ! sets value to depth average value
    ! val1
    val1_par(1) = val_avg_vp     ! left point values
    val1_par(2) = val_avg_vs
    val1_par(3) = val_avg_rho
    ! val2
    val2_par(1) = val_avg_vp     ! right point values
    val2_par(2) = val_avg_vs
    val2_par(3) = val_avg_rho
  else if (ix < 1 .and. ixplus <= Nx) then
    ! no point to the left, takes value from right one
    ! puts "virtual" left to border
    ix = 1
    ! sets value to value from right point
    ! val1
    if (FILL_MISSING_BORDER_POINTS_WITH_AVERAGE_VALUES) then
      val1_par(1) = val_avg_vp     ! left border
      val1_par(2) = val_avg_vs
      val1_par(3) = val_avg_rho
    else
      ! extend right point values
      val1_par(1) = EMC_vp(ixplus,ilat,idep)
      val1_par(2) = EMC_vs(ixplus,ilat,idep)
      val1_par(3) = EMC_rho(ixplus,ilat,idep)
    endif
    ! val2
    val2_par(1) = EMC_vp(ixplus,ilat,idep)
    val2_par(2) = EMC_vs(ixplus,ilat,idep)
    val2_par(3) = EMC_rho(ixplus,ilat,idep)
  else if (ix >= 1 .and. ixplus > Nx) then
    ! no point to the right, takes value from left one
    ! puts "virtual" right to border
    ixplus = Nx
    ! sets value to value from left point
    ! val1
    val1_par(1) = EMC_vp(ix,ilat,idep)
    val1_par(2) = EMC_vs(ix,ilat,idep)
    val1_par(3) = EMC_rho(ix,ilat,idep)
    ! val2
    if (FILL_MISSING_BORDER_POINTS_WITH_AVERAGE_VALUES) then
      val2_par(1) = val_avg_vp     ! right border
      val2_par(2) = val_avg_vs
      val2_par(3) = val_avg_rho
    else
      ! extend left point values
      val2_par(1) = EMC_vp(ix,ilat,idep)
      val2_par(2) = EMC_vs(ix,ilat,idep)
      val2_par(3) = EMC_rho(ix,ilat,idep)
    endif
  else
    ! left/right points valid
    ! sets value
    val1_par(1) = EMC_vp(ix,ilat,idep)
    val1_par(2) = EMC_vs(ix,ilat,idep)
    val1_par(3) = EMC_rho(ix,ilat,idep)

    val2_par(1) = EMC_vp(ixplus,ilat,idep)
    val2_par(2) = EMC_vs(ixplus,ilat,idep)
    val2_par(3) = EMC_rho(ixplus,ilat,idep)
  endif

  ! latitude
  if (iy < 1 .and. iyplus > Ny) then
    ! no valid point on this lat
    ! puts "virtual" points to border
    iy = 1
    iyplus = Ny
    ! sets value to depth average value
    ! val3
    val3_par(1) = val_avg_vp     ! lower point values
    val3_par(2) = val_avg_vs
    val3_par(3) = val_avg_rho
    ! val4
    val4_par(1) = val_avg_vp     ! upper point values
    val4_par(2) = val_avg_vs
    val4_par(3) = val_avg_rho
  else if (iy < 1 .and. iyplus <= Ny) then
    ! no point down, takes value from upper one
    ! puts "virtual" lower to border
    iy = 1
    ! sets value to value from right point
    ! val3
    if (FILL_MISSING_BORDER_POINTS_WITH_AVERAGE_VALUES) then
      val3_par(1) = val_avg_vp     ! right border
      val3_par(2) = val_avg_vs
      val3_par(3) = val_avg_rho
    else
      val3_par(1) = EMC_vp(ilon,iyplus,idep)
      val3_par(2) = EMC_vs(ilon,iyplus,idep)
      val3_par(3) = EMC_rho(ilon,iyplus,idep)
    endif
    ! val4
    val4_par(1) = EMC_vp(ilon,iyplus,idep)
    val4_par(2) = EMC_vs(ilon,iyplus,idep)
    val4_par(3) = EMC_rho(ilon,iyplus,idep)
  else if (iy >= 1 .and. iyplus > Ny) then
    ! no point up, takes value from lower one
    ! puts "virtual" upper to border
    iyplus = Ny
    ! sets value to value from left point
    ! val3
    val3_par(1) = EMC_vp(ilon,iy,idep)
    val3_par(2) = EMC_vs(ilon,iy,idep)
    val3_par(3) = EMC_rho(ilon,iy,idep)
    ! val4
    if (FILL_MISSING_BORDER_POINTS_WITH_AVERAGE_VALUES) then
      val4_par(1) = val_avg_vp     ! right border
      val4_par(2) = val_avg_vs
      val4_par(3) = val_avg_rho
    else
      val4_par(1) = EMC_vp(ilon,iy,idep)
      val4_par(2) = EMC_vs(ilon,iy,idep)
      val4_par(3) = EMC_rho(ilon,iy,idep)
    endif
  else
    ! lower/upper points valid
    ! sets value
    ! val3
    val3_par(1) = EMC_vp(ilon,iy,idep)
    val3_par(2) = EMC_vs(ilon,iy,idep)
    val3_par(3) = EMC_rho(ilon,iy,idep)
    ! val4
    val4_par(1) = EMC_vp(ilon,iyplus,idep)
    val4_par(2) = EMC_vs(ilon,iyplus,idep)
    val4_par(3) = EMC_rho(ilon,iyplus,idep)
  endif

  ! check index bounds
  if (ix < 1 .or. ix > Nx    ) stop 'Invalid interpolation point 1'
  if (ixplus < 1 .or. ixplus > Nx) stop 'Invalid interpolation point 2'
  if (iy < 1 .or. iy > Ny    ) stop 'Invalid interpolation point 3'
  if (iyplus < 1 .or. iyplus > Ny) stop 'Invalid interpolation point 4'

  if (ix == ixplus) stop 'Invalid interpolation point 1 and point 2'
  if (iy == iyplus) stop 'Invalid interpolation point 3 and point 4'

  ! interpolation
  dx = 1.d0
  dy = 1.d0
  gamma_interp_x = 0.d0
  gamma_interp_y = 0.d0

  ! grid point position
  lon = EMC_lon(ilon)
  lat = EMC_lat(ilat)

  ! determine ratio factor
  ! distance between the two interpolation points
  if (ix /= ixplus) dx = EMC_lon(ixplus) - EMC_lon(ix)
  ! ratio
  gamma_interp_x = (lon - dble(EMC_lon(ix))) / dx

  ! distance between the two interpolation points
  if (iy /= iyplus) dy = EMC_lat(iyplus) - EMC_lat(iy)
  ! ratio
  gamma_interp_y = (lat - dble(EMC_lat(iy))) / dy

  ! limits interpolation factors to be in range [0,1]
  if (gamma_interp_x < 0.d0) gamma_interp_x = 0.d0
  if (gamma_interp_x > 1.d0) gamma_interp_x = 1.d0
  if (gamma_interp_y < 0.d0) gamma_interp_y = 0.d0
  if (gamma_interp_y > 1.d0) gamma_interp_y = 1.d0

  ! bilinear interpolation
  vp_interp = interpolate_bilinear(val1_par(1),val2_par(1),val3_par(1),val4_par(1),gamma_interp_x,gamma_interp_y)
  vs_interp = interpolate_bilinear(val1_par(2),val2_par(2),val3_par(2),val4_par(2),gamma_interp_x,gamma_interp_y)
  rho_interp= interpolate_bilinear(val1_par(3),val2_par(3),val3_par(3),val4_par(3),gamma_interp_x,gamma_interp_y)

  end subroutine do_bilinear_interpolation

  !------------------------------------

  function interpolate_bilinear(val1,val2,val3,val4, &
                                gamma_interp_x,gamma_interp_y) result (interp_val)

  implicit none
  double precision, intent(in) :: val1,val2,val3,val4
  double precision, intent(in) :: gamma_interp_x,gamma_interp_y
  ! return value
  double precision :: interp_val

  ! interpolation rule
  ! uses bilinear interpolation
  interp_val =  &
       val1 * (1.d0-gamma_interp_x) * (1.d0-gamma_interp_y) + &
       val2 * gamma_interp_x        * (1.d0-gamma_interp_y) + &
       val3 * gamma_interp_x        * gamma_interp_y        + &
       val4 * (1.d0-gamma_interp_x) * gamma_interp_y

  end function interpolate_bilinear

  !------------------------------------

  subroutine do_shepard_interpolation(vp_interp,vs_interp,rho_interp)

! inverse distance weighted (IDL) interpolation (Shepard's method)

  use constants, only: HUGEVAL

  implicit none

  double precision, intent(inout) :: vp_interp,vs_interp,rho_interp

  ! local parameters
  double precision :: lon_ref,lat_ref,lon,lat,dlon,dlat
  double precision :: dist,dist_sq
  integer :: i,j,ix,iy
  integer :: icount
  ! inverse weighting
  double precision :: sum_weight,weight
  double precision :: sum_val_vp,sum_val_vs,sum_val_rho
  ! tapering
  double precision :: taper_val
  double precision, external :: cosine_taper

  ! reference grid point position
  lon_ref = EMC_lon(ilon)
  lat_ref = EMC_lat(ilat)

  ! initializes weighted values
  sum_val_vp = 0.d0
  sum_val_vs = 0.d0
  sum_val_rho = 0.d0

  sum_weight = 0.d0
  icount = 0

  ! finds defined points within a grid search
  do j = -NUMSTEPS_IY,NUMSTEPS_IY
    do i = -NUMSTEPS_IX,NUMSTEPS_IX
      ! grid index
      ix = ilon + i
      iy = ilat + j

      ! checks bounds
      if (ix < 1) cycle
      if (ix > Nx) cycle
      if (iy < 1) cycle
      if (iy > Ny) cycle

      ! skip missing points
      if (EMC_mask(ix,iy,idep)) cycle

      ! defined point location
      lon = EMC_lon(ix)
      lat = EMC_lat(iy)

      ! distance - as simple as it gets in degrees
      ! no Haversine function (that could be used for more accurate results, but more costly)
      dlon = lon_ref - lon
      dlat = lat_ref - lat
      dist_sq = dlon*dlon + dlat*dlat

      ! weight based on inverse distance
      weight = 1.d0 / dist_sq ! * taper_val

      ! Accumulate weighted sum
      sum_weight = sum_weight + weight

      ! computes taper
      dist = sqrt(dist_sq)
      taper_val = cosine_taper(dist, taper_dist_maximum)

      ! nearest point value w/ tapering
      vp_interp = EMC_vp(ix,iy,idep) * taper_val + (1.d0 - taper_val) * val_avg_vp
      vs_interp = EMC_vs(ix,iy,idep) * taper_val + (1.d0 - taper_val) * val_avg_vs
      rho_interp = EMC_rho(ix,iy,idep) * taper_val + (1.d0 - taper_val) * val_avg_rho

      sum_val_vp = sum_val_vp + weight * vp_interp
      sum_val_vs = sum_val_vs + weight * vs_interp
      sum_val_rho = sum_val_rho + weight * rho_interp

      ! counter
      icount = icount + 1
    enddo
  enddo

  !debug
  !print *,'debug: weighted interpolation: icount ',icount,' ilon/ilat ',ilon,ilat,'Nx/Ny',Nx,Ny,'sum_w',sum_weight
  !print *,'       reference lon/lat',lon_ref,lat_ref
  !print *,'       EMC lon min/max',EMC_lon_min,EMC_lon_max
  !print *,'       EMC lat min/max',EMC_lat_min,EMC_lat_max

  ! checks if we found a point
  if (icount == 0) then
    ! no near point found
    ! return background/average value
    vp_interp = val_avg_vp
    vs_interp = val_avg_vs
    rho_interp = val_avg_rho
  else
    ! Compute interpolated value
    vp_interp = sum_val_vp / sum_weight
    vs_interp = sum_val_vs / sum_weight
    rho_interp = sum_val_rho / sum_weight
  endif

  end subroutine do_shepard_interpolation

  end subroutine fill_EMC_missing_values_interpolated

!
!-------------------------------------------------------------------------------------------
!

  function cosine_taper(distance, max_distance) result(taper_val)

  use constants, only: TINYVAL,PI
  implicit none
  double precision,intent(in) :: distance, max_distance
  ! return result
  double precision :: taper_val

  ! initializes
  taper_val = 1.d0

  ! check
  if (abs(max_distance) < TINYVAL) return

  if (distance <= 0.d0) then
    taper_val = 1.d0
  else if (distance >= max_distance) then
    taper_val = 0.d0
  else
    taper_val = 0.5d0 * (1.d0 + cos(PI * distance / max_distance))
  endif

  end function cosine_taper


!-------------------------------------------------------------------------------------------
!
! main routines
!
!-------------------------------------------------------------------------------------------

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

  call bcast_all_singlei(EMC_depth_reference_level)

  ! broadcasts grid coordinates
  call bcast_all_singlei(EMC_lonlen)
  call bcast_all_singlei(EMC_latlen)
  call bcast_all_singlei(EMC_deplen)

  lonlen = EMC_lonlen
  latlen = EMC_latlen
  deplen = EMC_deplen

  ! allocate grid arrays for all other processes
  if (myrank /= 0) then
    allocate(EMC_lon(lonlen), &
             EMC_lat(latlen), &
             EMC_dep(deplen),stat=ier)
    if (ier /= 0) stop 'Error allocating lat,lon,dep arrays'
    EMC_lon(:) = 0.0; EMC_lat(:) = 0.0; EMC_dep(:) = 0.0
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
    ! model parameters
    allocate(EMC_vp(nx,ny,nz), &
             EMC_vs(nx,ny,nz), &
             EMC_rho(nx,ny,nz), stat=ier)
    if (ier /= 0) stop 'Error allocating EMC model arrays'
    EMC_vp(:,:,:) = 0.0; EMC_vs(:,:,:) = 0.0; EMC_rho(:,:,:) = 0.0

    ! mask for missing values
    allocate(EMC_mask(nx,ny,nz),stat=ier)
    if (ier /= 0) stop 'Error allocating mask'
    EMC_mask(:,:,:) = .false.

    ! average 1D model
    allocate(EMC_avg_rho(deplen),EMC_avg_vp(deplen),EMC_avg_vs(deplen),stat=ier)
    if (ier /= 0) stop 'Error allocating 1d average arrays'
    EMC_avg_rho(:) = 0.0; EMC_avg_vp(:) = 0.0; EMC_avg_vs(:) = 0.0
  endif

  ! broadcast the velocity model
  call bcast_all_cr(EMC_vp,nx * ny * nz)
  call bcast_all_cr(EMC_vs,nx * ny * nz)
  call bcast_all_cr(EMC_rho,nx * ny * nz)

  ! mask
  call bcast_all_l(EMC_mask,nx * ny * nz)

  ! average 1D model
  call bcast_all_cr(EMC_avg_vp,deplen)
  call bcast_all_cr(EMC_avg_vs,deplen)
  call bcast_all_cr(EMC_avg_rho,deplen)

  end subroutine model_emc_broadcast

!
!-------------------------------------------------------------------------------------------
!

  subroutine read_emc_model()

  use constants, only: IOUT

  use shared_parameters, only: REGIONAL_MESH_CUTOFF,REGIONAL_MESH_CUTOFF_DEPTH, &
    CENTER_LATITUDE_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH, &
    ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES, &
    NCHUNKS

  ! tapering
  use constants, only: TWO_PI,R_UNIT_SPHERE,DEGREES_TO_RADIANS,RADIANS_TO_DEGREES
  use shared_parameters, only: NCHUNKS,NEX_XI,NEX_ETA,ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES

  use model_emc_par

  implicit none

  ! local parameters
  integer :: ncid, status
  integer :: latid, lonid, depid    ! dimension ids
  integer :: latlen, lonlen, deplen ! dimension lengths
  integer :: dimlens(3)             ! dimension lengths (for array allocation purposes)
  integer :: varorderdims(3)        ! variable dimension order
  integer :: nx,ny,nz

  integer :: varid_vp, varid_vs, varid_rho    ! variable ids
  integer :: varid_lat, varid_lon, varid_dep
  integer :: dir_dep,dir
  integer :: i,ier

  real(kind=CUSTOM_REAL) :: missing_val_vp,missing_val_vs,missing_val_rho,missing_val_dep
  real(kind=CUSTOM_REAL) :: vp_min,vp_max,vs_min,vs_max,rho_min,rho_max
  real(kind=CUSTOM_REAL) :: dx,dx0,tmp_val
  character (len=MAX_STRING_LEN) :: filename

  ! chunk corners
  double precision :: corners_lat(4),corners_lon(4)

  ! tapering
  double precision :: ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD
  double precision :: element_size

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
    tmp_val = EMC_dep_max
    EMC_dep_max = - EMC_dep_min        ! switches sign, and also min/max
    EMC_dep_min = - tmp_val
    EMC_dep_dir = 2                    ! down
  endif

  ! flag for regional/global models
  EMC_is_regional = .true.
  if (abs((EMC_lat_max - EMC_lat_min) - 180.0) < TINYVAL .and. &
      abs((EMC_lon_max - EMC_lon_min) - 360.0) < TINYVAL) then
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
  if (FILL_MISSING_VALUES_WITH_AVERAGE) then
    write(IMAIN,*) '           using average 1D velocities for missing values'
  else if (FILL_MISSING_VALUES_WITH_INTERPOLATION) then
    write(IMAIN,*) '           using interpolated velocities for missing values'
  else
    write(IMAIN,*) '           using default background (PREM) values for missing values'
  endif
  if (EXTEND_EMC_MODEL_RANGE) then
    write(IMAIN,*) '           extending model range (laterally & vertically)'
  endif
  if (TAPER_OUTSIDE_RANGE) then
    write(IMAIN,*) '           using taper for outside points'
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

  ! velocity model
  if (.not. SCALE_MODEL) then
    ! no scaling
    ! must have complete vp,vs,rho setting
    if (varid_vp == 0) stop 'Error vp array variable not found'
    if (varid_vs == 0) stop 'Error vs array variable not found'
    if (varid_rho == 0) stop 'Error rho array variable not found'
  else
    ! allows for missing parameter scaling
    if (varid_vp /= 0 .and. (varid_vs == 0 .or. varid_rho == 0)) then
      if (varid_vs == 0) then
        write(IMAIN,*) '           scaling vs  from vp'
      endif
      if (varid_rho == 0) then
        write(IMAIN,*) '           scaling rho from vp'
      endif
      write(IMAIN,*)
      call flush_IMAIN()
    else if (varid_vs /= 0 .and. (varid_vp == 0 .or. varid_rho == 0)) then
      if (varid_vp == 0) then
        write(IMAIN,*) '           scaling vp  from vs'
      endif
      if (varid_rho == 0) then
        write(IMAIN,*) '           scaling rho from vp'
      endif
      write(IMAIN,*)
      call flush_IMAIN()
    else if (varid_vp /= 0 .and. varid_vs /= 0 .and. varid_rho /= 0) then
      write(IMAIN,*) '           model (vp,vs,rho) is complete'
      write(IMAIN,*)
      call flush_IMAIN()
    else
      write(IMAIN,*) '           model (vp,vs,rho) is incomplete:'
      if (varid_vp == 0)  write(IMAIN,*) '           vp is missing'
      if (varid_vs == 0)  write(IMAIN,*) '           vs is missing'
      if (varid_rho == 0) write(IMAIN,*) '           rho is missing'
      write(IMAIN,*)
      write(IMAIN,*) 'Please check your EMC model file "model.nc", and/or add more scaling relations to EMC models'
      call flush_IMAIN()
      stop 'Error velocity model is incomplete'
    endif
  endif

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

  ! checks up/down direction for depth
  call check_variable_attributes(ncid, varid_dep, EMC_dep_unit, dir_dep, missing_val_dep)

  ! converts depth to km
  ! units: 1==m, 2==km, 3==m/s, 4==km/s, 5==g/cm^3, 6==kg/cm^3, 7==kg/m^3
  if (EMC_dep_unit == 1) then
    ! converts to km
    EMC_dep(:) = EMC_dep(:) / 1000.0_CUSTOM_REAL
    EMC_dep_unit = 2    ! in km
  endif
  ! converts depth reference direction to positive being down (positive depth below sealevel, negative depth above)
  if (dir_dep == 1) then
    EMC_dep(:) = - EMC_dep(:)
    dir_dep = 2     ! down
  endif

  ! debug
  !print *,'debug: array lat min/max = ',minval(EMC_lat),maxval(EMC_lat),'header',EMC_lat_min,EMC_lat_max
  !print *,'debug: array lon min/max = ',minval(EMC_lon),maxval(EMC_lon),'header',EMC_lon_min,EMC_lon_max
  !print *,'debug: array dep min/max = ',minval(EMC_dep),maxval(EMC_dep),'header',EMC_dep_min,EMC_dep_max

  ! checks consistency between grid values and header infos
  if (abs(minval(EMC_lat) - EMC_lat_min) > 1.e-2) stop 'Error invalid EMC grid lat min'
  if (abs(maxval(EMC_lat) - EMC_lat_max) > 1.e-2) stop 'Error invalid EMC grid lat max'
  if (abs(minval(EMC_lon) - EMC_lon_min) > 1.e-2) stop 'Error invalid EMC grid lon min'
  if (abs(maxval(EMC_lon) - EMC_lon_max) > 1.e-2) stop 'Error invalid EMC grid lon max'
  if (abs(minval(EMC_dep) - EMC_dep_min) > 1.e-2) stop 'Error invalid EMC grid dep min'
  if (abs(maxval(EMC_dep) - EMC_dep_max) > 1.e-2) stop 'Error invalid EMC grid dep max'

  ! sets actual grid min/max values
  ! (header values might be truncated for readability as is the case in the SCEC CVM model file)
  EMC_lat_min = minval(EMC_lat)
  EMC_lat_max = maxval(EMC_lat)
  EMC_lon_min = minval(EMC_lon)
  EMC_lon_max = maxval(EMC_lon)
  EMC_dep_min = minval(EMC_dep)
  EMC_dep_max = maxval(EMC_dep)

  ! checks regular gridding
  ! lat & lon are usually regular, depth has often variable grid stepping (finer gridding in crust, coarser in mantle)
  ! lat
  EMC_regular_grid_lat = .true.
  EMC_dlat = (EMC_lat_max - EMC_lat_min) / (latlen-1)
  if (EMC_lat(2) < EMC_lat(1)) EMC_dlat = - EMC_dlat ! flips direction
  dx0 = EMC_lat(2) - EMC_lat(1)
  do i = 2,latlen-1
    dx = EMC_lat(i+1) - EMC_lat(i)
    if (abs(dx - dx0) > abs(1.e-2 * dx0) .and. abs(dx - EMC_dlat) > abs(1.e-2 * EMC_dlat)) then
      !debug
      if (VERBOSE) print *,'  irrregular gridding: dlat ',EMC_dlat,' - dx,dx0 = ',dx,dx0
      if (VERBOSE) print *
      ! set flag
      EMC_regular_grid_lat = .false.
      exit
    endif
    dx0 = dx
  enddo
  ! lon
  EMC_regular_grid_lon = .true.
  EMC_dlon = (EMC_lon_max - EMC_lon_min) / (lonlen-1)
  if (EMC_lon(2) < EMC_lon(1)) EMC_dlon = - EMC_dlon ! flips direction
  dx0 = EMC_lon(2) - EMC_lon(1)
  do i = 2,lonlen-1
    dx = EMC_lon(i+1) - EMC_lon(i)
    if (abs(dx - dx0) > abs(1.e-2 * dx0) .and. abs(dx - EMC_dlon) > abs(1.e-2 * EMC_dlon)) then
      !debug
      if (VERBOSE) print *,'  irrregular gridding: dlon ',EMC_dlon,' - dx,dx0 = ',dx,dx0
      if (VERBOSE) print *
      ! set flag
      EMC_regular_grid_lon = .false.
      exit
    endif
    dx0 = dx
  enddo
  ! depth
  EMC_regular_grid_dep = .true.
  EMC_ddep = (EMC_dep_max - EMC_dep_min) / (deplen-1)
  if (EMC_dep(2) < EMC_dep(1)) EMC_ddep = - EMC_ddep ! flips direction
  dx0 = EMC_dep(2) - EMC_dep(1)
  do i = 2,deplen-1
    dx = EMC_dep(i+1) - EMC_dep(i)
    if (abs(dx - dx0) > abs(1.e-2 * dx0) .and. abs(dx - EMC_ddep) > abs(1.e-2 * EMC_ddep)) then
      ! debug
      if (VERBOSE) print *,'  irrregular gridding: ddep ',EMC_ddep,' - dx,dx0 = ',dx,dx0
      if (VERBOSE) print *
      ! set flag
      EMC_regular_grid_dep = .false.
      exit
    endif
    dx0 = dx
  enddo

  ! grid spacing
  if (EMC_regular_grid_lat) then
    EMC_dlat = EMC_dlat
  else
    EMC_dlat = EMC_lat(2) - EMC_lat(1)
  endif
  if (EMC_regular_grid_lon) then
    EMC_dlon = EMC_dlon
  else
    EMC_dlon = EMC_lon(2) - EMC_lon(1)
  endif
  if (EMC_regular_grid_dep) then
    EMC_ddep = EMC_ddep
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
  if (EMC_regular_grid_lon) then
    write(IMAIN,*) '           increments dlon = ',EMC_dlon
  else
    write(IMAIN,*) '           increments dlon = variable (',EMC_dlon,',..)'
  endif
  if (EMC_regular_grid_lat) then
    write(IMAIN,*) '           increments dlat = ',EMC_dlat
  else
    write(IMAIN,*) '           increments dlat = variable (',EMC_dlat,',..)'
  endif
  if (EMC_regular_grid_dep) then
    write(IMAIN,*) '           increments ddep = ',EMC_ddep
  else
    write(IMAIN,*) '           increments ddep = variable (',EMC_ddep,',..)'
  endif
  write(IMAIN,*)
  if (EMC_depth_reference_level == 1) then
    write(IMAIN,*) '           depth reference level: earth surface'
  else if (EMC_depth_reference_level == 2) then
    write(IMAIN,*) '           depth reference level: sea level'
  else
    write(IMAIN,*) '           depth reference level: not recognized - assuming earth surface'
  endif
  write(IMAIN,*)
  call flush_IMAIN()

  ! model variables dimensions
  ! Check that the variable's dimensions are in the correct order
  if (varid_vp /= 0) then
    call check_dimorder(ncid, varid_vp, latid, lonid, depid, varorderdims)
  else if (varid_vs /= 0) then
    call check_dimorder(ncid, varid_vs, latid, lonid, depid, varorderdims)
  else if (varid_rho /= 0) then
    call check_dimorder(ncid, varid_rho, latid, lonid, depid, varorderdims)
  else
    print *,'Error: no parameter vp, vs or rho found for checking array dimensions'
    stop 'Error no model parameter vp, vs or rho found'
  endif

  ! user output
  write(IMAIN,*) '  order  : lon id ',lonid,' array index = ',varorderdims(lonid)
  write(IMAIN,*) '           lat id ',latid,' array index = ',varorderdims(latid)
  write(IMAIN,*) '           dep id ',depid,' array index = ',varorderdims(depid)
  write(IMAIN,*)
  call flush_IMAIN()

  ! safety check with stored ids
  if (latid == 0) stop 'Invalid lat ids between dimension and array infos'
  if (lonid == 0) stop 'Invalid lon ids between dimension and array infos'
  if (depid == 0) stop 'Invalid depth ids between dimension and array infos'

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
  write(IMAIN,*) '  size   : number of points lat         = ',latlen
  write(IMAIN,*) '           number of points lon         = ',lonlen
  write(IMAIN,*) '           number of points dep         = ',deplen
  write(IMAIN,*)
  write(IMAIN,*) '           array dimensions nx/ny/nz    = ',nx,ny,nz
  write(IMAIN,*) '           memory required per process  = ', &
                             sngl(3*dble(nx)*dble(ny)*dble(nz)*dble(CUSTOM_REAL)/1024./1024.),'MB'
  write(IMAIN,*)
  call flush_IMAIN()

  ! allocates model arrays
  allocate(EMC_vp(nx,ny,nz), &
           EMC_vs(nx,ny,nz), &
           EMC_rho(nx,ny,nz),stat=ier)
  if (ier /= 0) stop 'Error allocating vp,vs,rho arrays'
  EMC_vp(:,:,:) = 0.0; EMC_vs(:,:,:) = 0.0; EMC_rho(:,:,:) = 0.0

  ! Read vp, vs, rho
  if (varid_vp /= 0)  call check_status(nf90_get_var(ncid, varid_vp, EMC_vp))
  if (varid_vs /= 0)  call check_status(nf90_get_var(ncid, varid_vs, EMC_vs))
  if (varid_rho /= 0) call check_status(nf90_get_var(ncid, varid_rho, EMC_rho))

  ! gets units and missing values
  if (varid_vp /= 0) then
    call check_variable_attributes(ncid, varid_vp, EMC_vp_unit, dir, missing_val_vp)
    ! double-check as array might still have NaNs
    call convert_nan_to_missing(EMC_vp,missing_val_vp)
  endif
  if (varid_vs /= 0) then
    call check_variable_attributes(ncid, varid_vs, EMC_vs_unit, dir, missing_val_vs)
    ! double-check as array might still have NaNs
    call convert_nan_to_missing(EMC_vs,missing_val_vs)
  endif
  if (varid_rho /= 0) then
    call check_variable_attributes(ncid, varid_rho, EMC_rho_unit, dir, missing_val_rho)
    ! double-check as array might still have NaNs
    call convert_nan_to_missing(EMC_rho,missing_val_rho)
  endif

  ! Close netcdf file
  call check_status(nf90_close(ncid))

  ! scaling missing parameters
  if (SCALE_MODEL) then
    if (varid_vp /= 0 .and. (varid_vs == 0 .or. varid_rho == 0)) then
      ! VP provided
      ! Vs scaling
      if (varid_vs == 0) then
        ! user output
        write(IMAIN,*) '  scaling: Vs  from Vp   (Brocher scaling)'
        call flush_IMAIN()
        ! scales EMC_vs from EMC_vp
        call scale_Brocher_vs_from_vp()
        ! sets missing factor
        missing_val_vs = missing_val_vp
        where(abs(EMC_vp - missing_val_vp) < TINYVAL) EMC_vs = missing_val_vs
      endif
      ! Density scaling
      if (varid_rho == 0) then
        ! user output
        write(IMAIN,*) '  scaling: Rho from Vp   (Brocher scaling)'
        call flush_IMAIN()
        ! scales EMC_rho from EMC_vp
        call scale_Brocher_rho_from_vp()
        ! sets missing factor
        missing_val_rho = missing_val_vp
        where(abs(EMC_vp - missing_val_vp) < TINYVAL) EMC_rho = missing_val_rho
      endif
    else if (varid_vs /= 0 .and. (varid_vp == 0 .or. varid_rho == 0)) then
      ! VS provided
      ! Vp scaling
      if (varid_vp == 0) then
        ! user output
        write(IMAIN,*) '  scaling: Vp  from Vs   (Brocher scaling)'
        call flush_IMAIN()
        ! scales EMC_vs from EMC_vp
        call scale_Brocher_vp_from_vs()
        ! sets missing factor
        missing_val_vp = missing_val_vs
        where(abs(EMC_vs - missing_val_vs) < TINYVAL) EMC_vp = missing_val_vp
      endif
      ! Density scaling
      if (varid_rho == 0) then
        ! user output
        write(IMAIN,*) '  scaling: Rho from Vp   (Brocher scaling)'
        call flush_IMAIN()
        ! scales EMC_rho from EMC_vp
        call scale_Brocher_rho_from_vp()
        ! sets missing factor
        missing_val_rho = missing_val_vp
        where(abs(EMC_vp - missing_val_vp) < TINYVAL) EMC_rho = missing_val_rho
      endif
    else if (varid_vp /= 0 .and. varid_vs /= 0 .and. varid_rho /= 0) then
      ! complete, no scaling needed
      write(IMAIN,*) '  scaling: no scaling needed, all parameters provided'
      call flush_IMAIN()
      continue
    else
      stop 'Invalid model scaling relation not implemented yet'
    endif
    ! user output
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! mask missing values
  allocate(EMC_mask(nx,ny,nz),stat=ier)
  if (ier /= 0) stop 'Error allocating mask'
  EMC_mask(:,:,:) = .false.

  ! vp
  ! here mask is .false. for missing value, .true. for valid points
  !
  ! note: intel ifort compiler (2021.10.0 oneAPI) seems to have problems with the `where .. elsewhere ..` statement
  !       and crashes with a segmentation fault. however, it seems to work with only a `where ..` statement.
  !       as a work-around, we omit the `elsewhere` statement and initialize first the mask accordingly.
  !
  ! also, instead of comparing float values directly like `a == b`, we use an expression like `abs(a-b) < TINYVAL`
  ! to allow for some inaccuracy due to numerical precision.
  !
  !leads to ifort crashes:
  !where(EMC_vp == missing_val_vp)
  !  EMC_mask = .false.
  !elsewhere
  !  EMC_mask = .true.
  !end where
  !work-around:
  EMC_mask(:,:,:) = .true.
  where(abs(EMC_vp - missing_val_vp) < TINYVAL) EMC_mask = .false.
  ! min/max without missing values
  vp_min = minval(EMC_vp,mask=EMC_mask)
  vp_max = maxval(EMC_vp,mask=EMC_mask)

  ! vs
  ! here mask is .false. for missing value, .true. for valid points
  !leads to ifort crashes:
  !where(EMC_vs == missing_val_vs)
  !  EMC_mask = .false.
  !elsewhere
  !  EMC_mask = .true.
  !end where
  !work-around:
  EMC_mask(:,:,:) = .true.
  where(abs(EMC_vs - missing_val_vs) < TINYVAL) EMC_mask = .false.
  ! min/max without missing values
  vs_min = minval(EMC_vs,mask=EMC_mask)
  vs_max = maxval(EMC_vs,mask=EMC_mask)

  ! rho
  ! here mask is .false. for missing value, .true. for valid points
  !leads to ifort crashes:
  !where(EMC_rho == missing_val_rho)
  !  EMC_mask = .false.
  !elsewhere
  !  EMC_mask = .true.
  !end where
  !work-around:
  EMC_mask(:,:,:) = .true.
  where(abs(EMC_rho - missing_val_rho) < TINYVAL) EMC_mask = .false.
  ! min/max without missing values
  rho_min = minval(EMC_rho,mask=EMC_mask)
  rho_max = maxval(EMC_rho,mask=EMC_mask)

  ! overall mask
  ! opposite value as above: total mask is .true. for missing value to mask out those points
  !                          where at least one of the velocity model values (vp, vs or rho) is missing;
  !                          and mask==.false. for valid points.
  EMC_mask(:,:,:) = .false.
  where(abs(EMC_vp - missing_val_vp) < TINYVAL) EMC_mask = .true.
  where(abs(EMC_vs - missing_val_vs) < TINYVAL) EMC_mask = .true.
  where(abs(EMC_rho - missing_val_rho) < TINYVAL) EMC_mask = .true.

  ! user output
  write(IMAIN,*) '  model  : vp  min/max = ', vp_min,'/',vp_max
  write(IMAIN,*) '           vs  min/max = ', vs_min,'/',vs_max
  write(IMAIN,*) '           rho min/max = ', rho_min,'/',rho_max
  write(IMAIN,*)
  write(IMAIN,*) '           number of missing/incomplete model points is ',count(EMC_mask),' out of ',nx*ny*nz
  write(IMAIN,*) '                                                        ',(100.0*count(EMC_mask))/(nx*ny*nz),'%'
  write(IMAIN,*)
  call flush_IMAIN()

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
    write(IMAIN,*) '  simulation chunk:'
    write(IMAIN,*) '    center (lat,lon) at : (',sngl(CENTER_LATITUDE_IN_DEGREES),sngl(CENTER_LONGITUDE_IN_DEGREES),')'
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
        print *,'Error: Global simulation (NCHUNKS == 6) not supported for regional EMC models'
        stop 'EMC regional mesh incompatible with global simulations in Par_file'
      endif

      ! Par_file region must be within EMC grid region
      ! checks corners
      ! lat
      if (EMC_lat_min > minval(corners_lat)) then
        print *,'Error: Par_file defined region too large for EMC model.'
        print *,'Par_file : LATITUDE center = ',CENTER_LATITUDE_IN_DEGREES,' with chunk width = ', ANGULAR_WIDTH_ETA_IN_DEGREES
        print *,'           - > minimum lat = ',CENTER_LATITUDE_IN_DEGREES - 0.5 * ANGULAR_WIDTH_ETA_IN_DEGREES
        print *,'  chunck corners (lat,lon) at: (',corners_lat(1),corners_lon(1),'),(',corners_lat(2),corners_lon(2),')'
        print *,'                               (',corners_lat(3),corners_lon(3),'),(',corners_lat(4),corners_lon(4),')'
        print *,'EMC model: lat min        = ',EMC_lat_min
        print *,'Please choose a smaller region in Par_file for this EMC model.'
        stop 'Invalid region in Par_file for this EMC model.'
      endif
      if (EMC_lat_max < maxval(corners_lat)) then
        print *,'Error: Par_file defined region too large for EMC model.'
        print *,'Par_file : LATITUDE center = ',CENTER_LATITUDE_IN_DEGREES,' with chunk width = ', ANGULAR_WIDTH_ETA_IN_DEGREES
        print *,'           - > maximum lat = ',CENTER_LATITUDE_IN_DEGREES + 0.5 * ANGULAR_WIDTH_ETA_IN_DEGREES
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
        print *,'           - > minimum lon = ',CENTER_LONGITUDE_IN_DEGREES - 0.5 * ANGULAR_WIDTH_XI_IN_DEGREES
        print *,'  chunck corners (lat,lon) at: (',corners_lat(1),corners_lon(1),'),(',corners_lat(2),corners_lon(2),')'
        print *,'                               (',corners_lat(3),corners_lon(3),'),(',corners_lat(4),corners_lon(4),')'
        print *,'EMC model: lon min        = ',EMC_lon_min
        print *,'Please choose a smaller region in Par_file for this EMC model.'
        stop 'Invalid region in Par_file for this EMC model.'
      endif
      if (EMC_lon_max < maxval(corners_lon)) then
        print *,'Error: Par_file defined region too large for regional EMC models.'
        print *,'Par_file : LONGITUDE center = ',CENTER_LONGITUDE_IN_DEGREES,' with chunk width = ', ANGULAR_WIDTH_XI_IN_DEGREES
        print *,'           - > maximum lon = ',CENTER_LONGITUDE_IN_DEGREES + 0.5 * ANGULAR_WIDTH_XI_IN_DEGREES
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

  ! consistency conversions
  ! to have the model parameters always in the same unit (velocities in km/s, density in kg/m^3)

  ! converts density to default kg/m^3
  ! units: 1==m, 2==km, 3==m/s, 4==km/s, 5==g/cm^3, 6==kg/cm^3, 7==kg/m^3
  if (EMC_rho_unit == 5) then
    ! converts to kg/m^3
    ! rho [kg/m^3] = rho * 1000 [g/cm^3]
    EMC_rho(:,:,:) = EMC_rho(:,:,:) * 1000.0_CUSTOM_REAL
    EMC_rho_unit = 7 ! kg/m^3
  else if (EMC_rho_unit == 6) then
    ! converts to kg/m^3
    ! rho [kg/m^3] = rho * 1000000 [kg/cm^3]
    EMC_rho(:,:,:) = EMC_rho(:,:,:) * 1.e6_CUSTOM_REAL
    EMC_rho_unit = 7 ! kg/m^3
  endif
  ! converts velocity to default m/s
  if (EMC_vp_unit == 4) then
    ! converts to m/s
    EMC_vp(:,:,:) = EMC_vp(:,:,:) * 1000.0_CUSTOM_REAL
    EMC_vp_unit = 3
  endif
  if (EMC_vs_unit == 4) then
    ! converts to m/s
    EMC_vs(:,:,:) = EMC_vs(:,:,:) * 1000.0_CUSTOM_REAL
    EMC_vs_unit = 3
  endif

  ! re-orders arrays to have fixed lon/lat/depth indexing
  call reorder_EMC_model_arrays(lonid,latid,depid,varorderdims)

  ! model modifications
  !
  ! note: some EMC models might need some special modifications to avoid artificially low Vs values.
  !       to create an EMC model, some groups might interpolate velocities between fluid/air and solid grid points,
  !       which leads to very low Vs values, but still giving a solid medium, and not a fluid/air medium.
  !       we could try to mitigate these models a bit, but it's unclear what would be a reasonable modification.
  !
  !       for now, we take what we get and the user decides if that model is usable for seismic simulations...
  !
  !if (trim(EMC_model_id) == 'AlaskaFWANT_Vs2019') then
  !  ! AlaskaFWANT_Vs2019
  !  ! EMC file: YangAndGao-AlaskaFWANT-Vs2019-0.0.nc
  !  ! - the upper-most top Vs values look too artificial, either they are not valid (NaNf) or have
  !  !   very small, non-zero values (like 0.00313637 km/s, ..)
  !  !   -> we will replace this top layer with the values from the layer below
  !  !      however, this still doesn't solve very small Vs values as they progress even further down...
  !  ! user ouput
  !  write(IMAIN,*) '  using model modification for model: ',trim(EMC_model_id)
  !  write(IMAIN,*) '    replacing top layer with next layer values (to fix tiny Vs values)'
  !  write(IMAIN,*)
  !  call flush_IMAIN()
  !
  !  ! replaces first layer values with values from layer below
  !  !    top layer is at 0 km depth
  !  !    next         at 1.833 km
  !  !    next         at 3.666 km
  !  EMC_vp(:,:,1) = EMC_vp(:,:,2)
  !  EMC_vs(:,:,1) = EMC_vs(:,:,2)
  !  EMC_rho(:,:,1) = EMC_rho(:,:,2)
  !  EMC_mask(:,:,1) = EMC_mask(:,:,2)
  !  ! new model name
  !  EMC_model_id = trim(EMC_model_id) // '-modified'
  !endif

  ! extend solid surface values upwards
  call fill_EMC_surface_values()

  ! determine average 1D model
  call determine_EMC_depth_average()

  ! fills missing value grid points w/ 1d depth average values
  if (FILL_MISSING_VALUES_WITH_AVERAGE) then
    call fill_EMC_missing_values_average()
  endif

  ! fills missing value grid points w/ 1d depth average values
  if (FILL_MISSING_VALUES_WITH_INTERPOLATION) then
    call fill_EMC_missing_values_interpolated()
  endif

  ! tapering
  if (TAPER_OUTSIDE_RANGE) then
    ! compares fixed taper distances to mesh
    ! compute typical size of elements at the surface
    ! (normalized)
    if (NCHUNKS == 6) then
      ! estimation for global meshes (assuming 90-degree chunks)
      element_size = TWO_PI * R_UNIT_SPHERE / (4.d0 * NEX_XI)
    else
      ! estimation for 1-chunk meshes
      ANGULAR_WIDTH_XI_RAD = ANGULAR_WIDTH_XI_IN_DEGREES * DEGREES_TO_RADIANS
      ANGULAR_WIDTH_ETA_RAD = ANGULAR_WIDTH_ETA_IN_DEGREES * DEGREES_TO_RADIANS
      element_size = max( ANGULAR_WIDTH_XI_RAD/NEX_XI,ANGULAR_WIDTH_ETA_RAD/NEX_ETA ) * R_UNIT_SPHERE
    endif

    ! user output
    write(IMAIN,*) '  taper : horizontal taper distance (in deg) : ',taper_distance_maximum_hori_in_degrees,'(deg)'
    write(IMAIN,*) '          vertical   taper distance (in km)  : ',taper_distance_maximum_vert_in_km,'(km)'
    write(IMAIN,*)
    write(IMAIN,*) '          estimated element size (at surface)          : ',element_size * RADIANS_TO_DEGREES,'(deg)'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine read_emc_model

!
!-------------------------------------------------------------------------------------------
!

  subroutine model_EMC_crustmantle(iregion_code,r,theta,phi,vpv,vph,vsv,vsh,eta_aniso,rho)

  use constants
  use shared_parameters, only: R_PLANET,R_PLANET_KM,RHOAV,TOPOGRAPHY,ELLIPTICITY
  use meshfem_models_par, only: ibathy_topo
  use meshfem_models_par, only: rspl_ellip,ellipicity_spline,ellipicity_spline2,nspl_ellip

  use model_emc_par

  implicit none

  integer ,intent(in) :: iregion_code

  ! radius     - normalized by globe radius [0,1.x]
  ! theta/phi  - colatitude/longitude in rad (range theta/phi = [0,pi] / [0,2pi] (geocentric)
  double precision, intent(in) :: r,theta,phi

  ! absolute values, not perturbations
  double precision, intent(inout) :: vpv,vph,vsv,vsh,eta_aniso,rho

  ! local parameters
  double precision :: r0,lat,lon,colat,r_depth
  double precision :: vp_iso,vs_iso,rho_iso
  double precision :: vpl,vsl,rhol
  double precision :: scaleval_vel,scaleval_rho
  integer :: index_lat,index_lon,index_dep
  ! elevation
  double precision :: elevation
  ! flag for points with missing values
  logical :: is_point_outside

  ! tapering
  double precision :: dist_lon,dist_lat,dist,taper_val
  double precision, external :: cosine_taper

  ! checks if anything to do
  ! so far, EMC models only for crust/mantle are supported
  if (iregion_code /= IREGION_CRUST_MANTLE) return

  ! initializes local parameter values
  vpl = 0.d0
  vsl = 0.d0
  rhol = 0.d0

  ! corrects colatitude for ellipticity in case
  if (ELLIPTICITY) then
    ! converts to geographic colatitude
    call geocentric_2_geographic_colat_dble(theta,colat)
  else
    ! spherical mesh, geocentric and geographic colatitudes are the same
    colat = theta
  endif

  ! geographic lat/lon
  ! latitude in degrees in [-90,90]
  lat = (PI_OVER_TWO - colat) * RADIANS_TO_DEGREES
  ! longitude in degrees in [0,360]
  lon = phi * RADIANS_TO_DEGREES

  ! determines radius of earth surface at position lat/lon (to calculate depth of position r)
  ! normalized surface radius
  r0 = R_UNIT_SPHERE

  ! finds elevation of position
  if (TOPOGRAPHY) then
    ! EMC depth reference level 1==earth surface, 2==sea level),
    if (EMC_depth_reference_level == 1) then
      ! depth relative to earth surface
      ! gets elevation
      ! needs lat in range [-90,90], lon in range [0,360]
      call get_topo_bathy(lat,lon,elevation,ibathy_topo)
      ! adapts surface radius
      r0 = r0 + elevation/R_PLANET
    else
      ! depth relative to sea level
      ! with sea level supposed to be normalized at 1
      r0 = r0       ! already done
    endif
  endif

  ! ellipticity
  if (ELLIPTICITY) then
    ! adds ellipticity factor to radius (needs geocentric colatitude)
    call add_ellipticity_rtheta(r0,theta,nspl_ellip,rspl_ellip,ellipicity_spline,ellipicity_spline2)
  endif

  ! depth with respect to surface (in km)
  r_depth = R_PLANET_KM * (r0 - r)

  ! checks depth limits
  ! depth reference level 1==earth surface, 2==sea level),
  if (EMC_depth_reference_level == 1) then
    ! depth relative to earth surface
    ! checks limit (mesh point cannot be in air)
    if (r_depth < 0.d0) r_depth = 0.d0
  else
    ! depth relative to sea level
    ! checks limit (mesh point cannot be higher than EMC limit)
    ! depth positive towards earth center
    !if (r_depth < EMC_dep_min) r_depth = EMC_dep_min
    ! will be checked in any case below
    continue
  endif

  ! general depth check (mesh point cannot be higher than EMC limit)
  !
  ! note: this has also the effect that the EMC model velocities extend upwards, i.e., for EMC models that have been defined
  !       only below a given depth, the velocities from the upper most depth will be extended up to the surface.
  !       for example, the Ward (2018) model starts at a depth of 1 km down to a depth of 70 km. the velocities at 1 km depth
  !       will be used for any shallower depths.
  !
  ! depth positive towards earth center
  if (r_depth < EMC_dep_min) r_depth = EMC_dep_min

  ! note: here, the choice is that any point outside of the EMC model box falls back to its background velocity model value,
  !       which by default will be the isotropic PREM.
  !
  !       still, this could lead to a problem with the depth in case there is a discrepancy between the elevation of the point
  !       when topography is used and if the EMC model starts at zero depth.
  !
  !       for different models, depth refers either to depth with respect to sea level or surface topography?
  !       so, if depth refers to sea level, then the depth value is more like a z-coordinate value
  !       (e.g., -2.0 for 2km elevation, 10.0 to 10km depth below sea level).
  !       on the other hand, if depth refers to the surface, then the depth value is more like a relative value
  !       (e.g., 2.0 for 2km depth, 10.0 to 10km depth below surface).
  !
  !       thus, for the mesher, if we determine the point values based on the grid points of the initial, spherical grid,
  !       and then stretch the points to match topography, the point location search here includes no topography yet,
  !       and depth should always be relative to surface.
  !
  !       however, if depth given is with respect to sea level, we should start the point search only with the
  !       final point location after topography stretching.
  !       -> there would be a small discrepancy if we transform depth to be always relative to surface, read in the model
  !          based on the initial spherical mesh and then stretch between bottom to top to accommodate surface topography.
  !          the transformation of (elevation) depth to depth with respect to surface would need to account for the
  !          stretching factor.
  !
  !       todo: determine depth for initial or final grid points? left for future improvements...

  ! shift lon [0,360] if EMC longitudes start at negative lon
  if (EMC_lon_min < -180.0) then
    if (lon > 0.0) lon = lon - 360.0d0      ! range [-360,0]
  else if (EMC_lon_min < 0.0 .and. EMC_lon_max <= 0.0) then
    if (lon > 0.0) lon = lon - 360.0d0      ! range [-360,0]
  else if (EMC_lon_min < 0.0) then
    if (lon > 180.0) lon = lon - 360.0d0    ! range [-180,180]
  endif

  ! safety check
  ! we assume from here on that the increments EMC_dlon, EMC_dlat and EMC_ddep are strictly positive.
  ! this should have been done by re-orderering the arrays in the read routine if needed.
  if (EMC_dlon <= 0.0_CUSTOM_REAL) stop 'Invalid EMC dlon increment'
  if (EMC_dlat <= 0.0_CUSTOM_REAL) stop 'Invalid EMC dlat increment'
  if (EMC_ddep <= 0.0_CUSTOM_REAL) stop 'Invalid EMC ddep increment'

  ! finds grid index
  ! lat [-90,90]
  if (EMC_regular_grid_lat) then
    ! indexing starts at 1,..
    index_lat = floor((lat - EMC_lat(1))/EMC_dlat) + 1
  else
    ! gets index of value in EMC_lat closest to given lat
    index_lat = minloc(abs(EMC_lat(:) - lat),dim=1)
    ! determine closest, but smaller value (like floor(..) above) than given lat
    if (lat < EMC_lat(index_lat)) index_lat = index_lat - 1
  endif
  ! checks bounds to allow index_lat go beyond the range [1,latlen], indicating that the point lies outside
  if (lat < EMC_lat(1)) index_lat = 0
  if (lat > EMC_lat(EMC_latlen)) index_lat = EMC_latlen + 1

  ! lon
  if (EMC_regular_grid_lon) then
    ! indexing starts at 1,..
    index_lon = floor((lon - EMC_lon(1))/EMC_dlon) + 1
  else
    ! gets index of value in EMC_lon closest to given lon
    index_lon = minloc(abs(EMC_lon(:) - lon),dim=1)
    ! determine closest, but smaller value (like floor(..) above) than given lon
    if (lon < EMC_lon(index_lon)) index_lon = index_lon - 1
  endif
  ! checks bounds to allow index_lat go beyond the range [1,lonlen], indicating that the point lies outside
  if (lon < EMC_lon(1)) index_lon = 0
  if (lon > EMC_lon(EMC_lonlen)) index_lon = EMC_lonlen + 1

  ! depth
  if (EMC_regular_grid_dep) then
    ! indexing starts at 1,..
    index_dep = floor((r_depth - EMC_dep(1))/EMC_ddep) + 1
  else
    ! gets index of value in EMC_dep closest to given dep
    index_dep = minloc(abs(EMC_dep(:) - r_depth),dim=1)
    ! determine closest, but smaller value (like floor(..) above) than given depth
    ! for example: r_depth = 30.8 -> EMC_dep == 31 closest for minloc index
    !              we want index with EMC_dep == 30 and use the cell between [30,31]km depth for interpolation
    ! in case depth increases in array: ddep > 0 and index with smaller depth is -1
    ! in case depth decreases in array: ddep < 0 and index with smaller depth is +1
    if (r_depth < EMC_dep(index_dep)) index_dep = index_dep - 1
  endif
  ! checks bounds to allow index_dep go beyond the range [1,deplen], indicating that the point lies outside
  if (r_depth < EMC_dep(1)) index_dep = 0
  if (r_depth > EMC_dep(EMC_deplen)) index_dep = EMC_deplen + 1

  ! extend model laterally to points outside EMC region
  if (EXTEND_EMC_MODEL_RANGE) then
    ! limits index
    ! lat
    if (index_lat < 1) index_lat = 1                        ! interpolation will take care
    if (index_lat > EMC_latlen) index_lat = EMC_latlen
    ! lon
    if (index_lon < 1) index_lon = 1                        ! interpolation will take care
    if (index_lon > EMC_lonlen) index_lon = EMC_lonlen
    ! depth extensions
    ! note: we will always extend the model values upwards to the surface if needed;
    !       this will allow for variations on top.
    !       for deeper parts, extending the model smears out structures vertically (looks like plume-like structures below).
    !       one could also replace value below the model range with average/background velocities.
    !       however, this introduces an artificial discontinuity at the bottom of the EMC model.
    ! here, we choose to extend in both directions as we will use by default a taper to taper off values outside the model.
    if (index_dep < 1) index_dep = 1                        ! extends model to top
    if (index_dep > EMC_deplen) index_dep = EMC_deplen      ! extends model below
    ! or
    !!if (index_dep > EMC_deplen) index_dep = EMC_deplen +1   ! replaces bottom by average/background
  endif

  ! checks if point outside
  is_point_outside = .false.
  if (index_lat < 1 .or. index_lat > EMC_latlen) is_point_outside = .true.
  if (index_lon < 1 .or. index_lon > EMC_lonlen) is_point_outside = .true.
  if (index_dep > EMC_deplen) is_point_outside = .true.

  ! background velocities (from PREM)
  ! Voigt average of input (for masked points where model values are missing)
  vp_iso = sqrt( (2.d0*vpv*vpv + vph*vph)/3.d0 )
  vs_iso = sqrt( (2.d0*vsv*vsv + vsh*vsh)/3.d0 )
  rho_iso = rho

  ! input background values are non-dimensionalized already
  ! - here we put dimension back to m/s (and kg/m^3) to mix with EMC model values for interpolation
  vp_iso = vp_iso * R_PLANET * sqrt(PI*GRAV*RHOAV)  ! to m/s
  vs_iso = vs_iso * R_PLANET * sqrt(PI*GRAV*RHOAV)  ! to m/s
  rho_iso = rho_iso * RHOAV                         ! to kg/m3

  ! using average 1D model depth values
  if (FILL_MISSING_VALUES_WITH_AVERAGE .or. FILL_MISSING_VALUES_WITH_INTERPOLATION) then
    if (index_dep >= 1 .and. index_dep <= EMC_deplen) then
      ! takes average values at depth
      call get_interpolated_1d_depth_average(r_depth,index_dep,vp_iso,vs_iso,rho_iso)
    else
      ! point outside model range
      ! last layer values
      vp_iso = EMC_avg_vp(EMC_deplen)
      vs_iso = EMC_avg_vs(EMC_deplen)
      rho_iso = EMC_avg_rho(EMC_deplen)
    endif
  endif

  ! interpolation
  if (.not. is_point_outside) then
    ! takes model values interpolated from EMC model
    call get_interpolated_EMC_values(index_lon,index_lat,index_dep,lon,lat,r_depth,rho_iso,vp_iso,vs_iso,rhol,vpl,vsl)
  else
    ! set to average/background value
    rhol = rho_iso
    vpl = vp_iso
    vsl = vs_iso
  endif

  ! point outside model
  if (TAPER_OUTSIDE_RANGE) then
    ! determines if point outside based on lat/lon position
    ! get maximum distance in radians
    ! lon
    if (lon < EMC_lon(1)) then
      dist_lon = abs(EMC_lon(1) - lon)
    else if (lon > EMC_lon(EMC_lonlen)) then
      dist_lon = abs(lon - EMC_lon(EMC_lonlen))
    else
      dist_lon = 0.d0
    endif
    ! lat
    if (lat < EMC_lat(1)) then
      dist_lat = abs(EMC_lat(1) - lat)
    else if (lat > EMC_lat(EMC_latlen)) then
      dist_lat = abs(lat - EMC_lat(EMC_latlen))
    else
      dist_lat = 0.d0
    endif

    ! distance in degrees
    dist = max(dist_lon,dist_lat)

    ! apply horizontal taper
    if (dist > 0.d0) then
      ! position outside EMC model range
      ! gets cosine taper value
      taper_val = cosine_taper(dist,taper_distance_maximum_hori_in_degrees)

      ! taper between (interpolated) and background/average value
      vpl = vpl * taper_val + (1.d0 - taper_val) * vp_iso
      vsl = vsl * taper_val + (1.d0 - taper_val) * vs_iso
      rhol = rhol * taper_val + (1.d0 - taper_val) * rho_iso
    endif

    ! vertical distance
    if (r_depth > EMC_dep(EMC_deplen)) then
      dist = abs(r_depth - EMC_dep(EMC_deplen))
    else
      dist = 0.d0
    endif

    ! apply vertical taper
    if (dist > 0.d0) then
      ! position below EMC model range
      ! gets cosine taper value
      taper_val = cosine_taper(dist,taper_distance_maximum_vert_in_km)

      ! taper between (interpolated) and background/average value
      vpl = vpl * taper_val + (1.d0 - taper_val) * vp_iso
      vsl = vsl * taper_val + (1.d0 - taper_val) * vs_iso
      rhol = rhol * taper_val + (1.d0 - taper_val) * rho_iso
    endif
  endif

  !debug
  !if (r_depth > 10.1 .and. r_depth < 15.1) &
  !if (vsl <= 1.0) &
  !  print *,'debug: lat/lon/dep = ',lat,lon,r_depth,'vp/vs/rho = ',vpl,vsl,rhol,'iso vp/vs/rho',vp_iso,vs_iso,rho_iso, &
  !          'index lat/lon/dep',index_lat,index_lon,index_dep

  ! only uses solid domain values
  ! (no fluid domain such as oceans are modelled so far, effect gets approximated by ocean load)
  ! returns model values if non-zero
  if (vsl > TINYVAL) then
    ! non-dimensionalize
    scaleval_rho = 1.0d0 / RHOAV                             ! from kg/m3
    scaleval_vel = 1.0d0 / (R_PLANET * sqrt(PI*GRAV*RHOAV))  ! from m/s (scaleval_vel == 1.459769779014117E-004)

    rhol = rhol * scaleval_rho
    vpl = vpl * scaleval_vel
    vsl = vsl * scaleval_vel

    ! converts isotropic values to transverse isotropic
    vpv = vpl
    vph = vpl
    vsv = vsl
    vsh = vsl
    eta_aniso = 1.d0
    rho = rhol
  endif

  end subroutine model_EMC_crustmantle

!
!-------------------------------------------------------------------------------------------
!

  subroutine get_interpolated_1d_depth_average(r_depth,index_dep,vp_iso,vs_iso,rho_iso)

  use model_emc_par, only: EMC_dep,EMC_deplen, &
                           EMC_avg_vp,EMC_avg_vs,EMC_avg_rho

  implicit none

  double precision,intent(in) :: r_depth
  integer, intent(in) :: index_dep
  double precision, intent(inout) :: vp_iso,vs_iso,rho_iso

  ! local parameters
  ! 1d average
  double precision :: dz,interp_z

  if (index_dep < EMC_deplen) then
    ! interpolated between two depth values
    dz = EMC_dep(index_dep+1) - EMC_dep(index_dep)
    interp_z = (r_depth - dble(EMC_dep(index_dep))) / dz

    if (interp_z < 0.d0) interp_z = 0.d0
    if (interp_z > 1.d0) interp_z = 1.d0

    ! linear interpolation
    vp_iso  = EMC_avg_vp(index_dep) * (1.d0 - interp_z)  + EMC_avg_vp(index_dep+1) * interp_z
    vs_iso  = EMC_avg_vs(index_dep) * (1.d0 - interp_z)  + EMC_avg_vs(index_dep+1) * interp_z
    rho_iso = EMC_avg_rho(index_dep) * (1.d0 - interp_z) + EMC_avg_rho(index_dep+1) * interp_z
  else
    ! last layer values
    vp_iso = EMC_avg_vp(index_dep)
    vs_iso = EMC_avg_vs(index_dep)
    rho_iso = EMC_avg_rho(index_dep)
  endif

  end subroutine get_interpolated_1d_depth_average

!
!-------------------------------------------------------------------------------------------
!

  subroutine get_interpolated_EMC_values(index_lon,index_lat,index_dep,lon,lat,r_depth,rho_iso,vp_iso,vs_iso,rhol,vpl,vsl)

  use model_emc_par, only: EMC_rho,EMC_vp,EMC_vs, &
                           EMC_lat,EMC_lon,EMC_dep, &
                           EMC_lonlen,EMC_latlen,EMC_deplen

  implicit none

  integer, intent(in) :: index_lon,index_lat,index_dep
  double precision, intent(in) :: r_depth,lat,lon
  double precision, intent(in) :: vp_iso,vs_iso,rho_iso
  double precision, intent(inout) :: vpl,vsl,rhol

  ! local parameters
  integer :: ix,iy,iz
  integer :: Nx,Ny,Nz
  ! positioning/interpolation
  double precision :: gamma_interp_x,gamma_interp_y,gamma_interp_z
  double precision :: dx,dy,dz

  ! initializes
  dx = 1.d0
  dy = 1.d0
  dz = 1.d0
  gamma_interp_x = 0.d0
  gamma_interp_y = 0.d0
  gamma_interp_z = 0.d0

  ! sets array indices
  ! order: lon/lat/dep (re-ordered EMC model arrays should have this order)
  ix = index_lon
  iy = index_lat
  iz = index_dep

  Nx = EMC_lonlen
  Ny = EMC_latlen
  Nz = EMC_deplen

  if (ix >= 1 .and. ix <= Nx) then
    if (ix < Nx) dx = EMC_lon(ix+1) - EMC_lon(ix)
    gamma_interp_x = (lon - dble(EMC_lon(ix))) / dx
  endif
  if (iy >= 1 .and. iy <= Ny) then
    if (iy < Ny) dy = EMC_lat(iy+1) - EMC_lat(iy)
    gamma_interp_y = (lat - dble(EMC_lat(iy))) / dy
  endif
  if (iz >= 1 .and. iz <= Nz) then
    if (iz < Nz) dz = EMC_dep(iz+1) - EMC_dep(iz)
    gamma_interp_z = (r_depth - dble(EMC_dep(iz))) / dz
  endif

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

  !debug: checks gamma factors
  !if (gamma_interp_x < 0.d0 .or. gamma_interp_y < 0.d0 .or. gamma_interp_z < 0.d0 .or. &
  !    gamma_interp_x > 1.d0 .or. gamma_interp_y > 1.d0 .or. gamma_interp_z > 1.d0) then
  !  print *,'Error: position has invalid gamma interpolation factors: '
  !  print *,'  point index            : ',ix,iy,iz,'Nx/Ny/Nz',Nx,Ny,Nz
  !  print *,'  location (lon/lat/dep) : ',sngl(lon),sngl(lat),sngl(r_depth)
  !  print *,'  grid point             : ',EMC_lon(ix),EMC_lat(iy),EMC_dep(iz)
  !  print *,'  interpolation gamma    : ',gamma_interp_x,gamma_interp_y,gamma_interp_z
  !  print *,'  interpolation dx/dy/dz : ',dx,dy,dz
  !  !call exit_MPI(myrank,'Error point index in EMC model routine')
  !endif

  ! limits interpolation factors to be in range [0,1]
  if (gamma_interp_x < 0.d0) gamma_interp_x = 0.d0
  if (gamma_interp_x > 1.d0) gamma_interp_x = 1.d0
  if (gamma_interp_y < 0.d0) gamma_interp_y = 0.d0
  if (gamma_interp_y > 1.d0) gamma_interp_y = 1.d0
  if (gamma_interp_z < 0.d0) gamma_interp_z = 0.d0
  if (gamma_interp_z > 1.d0) gamma_interp_z = 1.d0

  ! model vp
  ! stores interpolated value
  vpl = interpolate_trilinear(EMC_vp,vp_iso,ix,iy,iz,gamma_interp_x,gamma_interp_y,gamma_interp_z)

  ! model vs
  ! stores interpolated value
  vsl = interpolate_trilinear(EMC_vs,vs_iso,ix,iy,iz,gamma_interp_x,gamma_interp_y,gamma_interp_z)

  ! model rho
  ! stores interpolated value
  rhol = interpolate_trilinear(EMC_rho,rho_iso,ix,iy,iz,gamma_interp_x,gamma_interp_y,gamma_interp_z)

contains

  function interpolate_trilinear(EMC_par,EMC_par_iso,ix,iy,iz, &
                                 gamma_interp_x,gamma_interp_y,gamma_interp_z) result (interp_val)

  use constants, only: CUSTOM_REAL

  use model_emc_par, only: EMC_mask

  implicit none
  real(kind=CUSTOM_REAL), dimension(Nx,Ny,Nz), intent(in) :: EMC_par
  double precision, intent(in) :: EMC_par_iso
  integer, intent(in) :: ix,iy,iz
  double precision, intent(in) :: gamma_interp_x,gamma_interp_y,gamma_interp_z
  ! return value
  double precision :: interp_val

  ! local variables
  double precision :: val1,val2,val3,val4,val5,val6,val7,val8

  ! model parameter (vp, vs or rho)
  val1 = EMC_par(ix  ,iy  ,iz  )
  val2 = EMC_par(ix+1,iy  ,iz  )
  val3 = EMC_par(ix+1,iy+1,iz  )
  val4 = EMC_par(ix  ,iy+1,iz  )
  val5 = EMC_par(ix  ,iy  ,iz+1)
  val6 = EMC_par(ix+1,iy  ,iz+1)
  val7 = EMC_par(ix+1,iy+1,iz+1)
  val8 = EMC_par(ix  ,iy+1,iz+1)

  if (EMC_mask(ix  ,iy  ,iz  )) val1 = EMC_par_iso  ! takes input parameter (coming from background model PREM)
  if (EMC_mask(ix+1,iy  ,iz  )) val2 = EMC_par_iso
  if (EMC_mask(ix+1,iy+1,iz  )) val3 = EMC_par_iso
  if (EMC_mask(ix  ,iy+1,iz  )) val4 = EMC_par_iso
  if (EMC_mask(ix  ,iy  ,iz+1)) val5 = EMC_par_iso
  if (EMC_mask(ix+1,iy  ,iz+1)) val6 = EMC_par_iso
  if (EMC_mask(ix+1,iy+1,iz+1)) val7 = EMC_par_iso
  if (EMC_mask(ix  ,iy+1,iz+1)) val8 = EMC_par_iso

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

  end function interpolate_trilinear

  end subroutine get_interpolated_EMC_values
