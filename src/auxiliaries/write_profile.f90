!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            August 2013
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

! This is a utility code written to ouput profiles through the model used by SPECFEM,
! as it were if it was sampled on a radial line through the Earth at 2x2 degree
! intervals, with discontinuities honored by the mesh and those of the MOHO and
! crust/ocean indicated by repeating points with different values.
!
! The code shortcuts through meshfem3D -> create_regions_mesh ->
! create_regular_elements -> compute_element_properties -> compute_element_properties ->
! get_model, cutting and pasting the relevant parts.
!
! I have tried to indicate where different parts come from to facilitate updates.
! - vala hjorleifsdottir (vala@geofisica.unam.mx)

! note:
! This program writes profiles of material properties for the model chosen in Par_file.
! It samples directly the original choosen input MODEL. It does not sample the GLL-mesh which is output by the mesher.
! Thus, there might be slight differences, e.g., in how TISO properties are used in the solver.
!
! Most of the code is from get_model.f90 and is marked with comments.
! Copy paste the corresponding part from the most recent get_model.f90 and paste it into write_profile.f90.
! Updated calls to subroutines may be necessary, in particular to read_compute_parameters which often changes.
!
! The code should write profiles by default every 2x2 degrees for the entire globe.
!
!July 15th, 2008, updated Oct 27th, 2010
!originally by: Vala Hjorleifsdottir
!
! usage:
! - default for 2x2 degree profiles around the globe
!   ./bin/xwrite_profile
!
! - for a 5x5 degree
!   ./bin/xwrite_profile 5 5
!
! - for a single profile at location lat/lon = 47.37/-25.0
!   ./bin/xwrite_profile 1 1 47.37 -25.0
!


  program xwrite_profile

  use constants, only: &
    GRAV,PI,RADIANS_TO_DEGREES,MAX_STRING_LEN,IMAIN,ISTANDARD_OUTPUT, &
    IREGION_INNER_CORE,IREGION_OUTER_CORE,IREGION_CRUST_MANTLE, &
    IFLAG_IN_FICTITIOUS_CUBE,IFLAG_INNER_CORE_NORMAL, &
    IFLAG_MIDDLE_CENTRAL_CUBE,IFLAG_TOP_CENTRAL_CUBE,IFLAG_BOTTOM_CENTRAL_CUBE, &
    IFLAG_OUTER_CORE_NORMAL, &
    IFLAG_MANTLE_NORMAL,IFLAG_670_220,IFLAG_220_80,IFLAG_80_MOHO,IFLAG_CRUST, &
    IPLANET_EARTH, &
    myrank

  use meshfem_models_par, only: &
    CRUSTAL,ONE_CRUST

  ! ellipticity
  use meshfem_models_par, only: nspl_ellip,rspl_ellip,ellipicity_spline,ellipicity_spline2

  use shared_parameters, only: &
    doubling_index,MAX_NUMBER_OF_MESH_LAYERS,rmaxs,rmins, &
    MODEL,MODEL_GLL,OUTPUT_FILES, &
    OCEANS,TOPOGRAPHY,ELLIPTICITY,TRANSVERSE_ISOTROPY, &
    PLANET_TYPE,R_PLANET,RCMB,RICB,RMOHO,RMOHO_FICTITIOUS_IN_MESHER,ROCEAN,RHOAV

  implicit none

!-------------------------------------------------------
! USER parameters

  ! initial position
  double precision,parameter :: COLAT_0 = 1.d0
  double precision,parameter :: LON_0   = 1.d0

  ! colatitude loop range (in degrees)
  integer,parameter :: COLAT_istart = 0  ! 0
  integer,parameter :: COLAT_iend   = 89 ! 89  for 2x2 increments - > [0,178]

  ! longitude loop range (in degrees)
  integer,parameter :: LON_istart   = 0   ! 0
  integer,parameter :: LON_iend     = 179 ! 179  for 2x2 increments - > [0,358]

  ! sampling intervals
  double precision,parameter :: INCR_lat = 2.d0 ! 2 degrees
  double precision,parameter :: INCR_lon = 2.d0 ! 2 degrees

  ! depth sampling
  ! every 1km sampling depth for layers above this depth, otherwise uses 10 km depth samples
  ! (note: for 3D models, the moho-80km layer gets stretched down to 120km)
  double precision,parameter :: DELTA_HIRES_DEPTH = 150.0d0 ! in km

!-------------------------------------------------------

  ! local parameters
  ! parameters needed to store the radii of the grid points in the spherically symmetric Earth
  integer :: idoubling
  integer :: NUMBER_OF_MESH_LAYERS
  ! for loop on all the slices
  integer :: iregion_code
  ! proc numbers for MPI
  integer :: sizeprocs

  ! arguments
  integer :: count
  character(len=MAX_STRING_LEN) :: arg
  logical :: use_userdefined_lat

  double precision :: rmin,rmax,rmax_last,r_prem
  double precision :: r,theta,phi
  double precision :: theta_degrees,phi_degrees
  double precision :: lat,lon,elevation
  double precision :: moho

  double precision :: rho,Qkappa,Qmu
  double precision :: vpv,vph,vsv,vsh,eta_aniso
  double precision :: scaleval

  integer :: iline,iline_icb,iline_cmb,iline_moho,iline_ocean
  integer :: nit,ilayer,idep
  integer :: i,j,ier

  integer :: icolat_start,icolat_end,ilon_start,ilon_end
  double precision :: initial_colat,initial_lon
  ! lat/lon increments
  double precision :: delta_lat,delta_lon
  ! depth increment
  double precision :: delta

  ! tolerance in checks for equal (==) float values
  double precision, parameter :: TOL_ZERO = 1.d-12

  character(len=MAX_STRING_LEN) :: outfile
  character(len=7) :: str_info

! ************** PROGRAM STARTS HERE **************

  ! initializes MPI to be able to call meshfem routines below
  ! starts mpi
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  ! checks number of processes
  ! note: must run as a single process, if necessary with: mpirun -np 1 ./bin/xwrite_profile
  if (sizeprocs /= 1) then
    ! usage info
    if (myrank == 0) then
      print *,'This program executes as sequential program.'
      print *,'Invalid number of processes used: ', sizeprocs, ' procs'
      print *
    endif
    call usage()
    call abort_mpi()
  endif

  ! main output file
  ! used by some of the mesher subroutines, only written to by process 0
  if (myrank == 0) then
    if (IMAIN /= ISTANDARD_OUTPUT) &
      open(unit=IMAIN,file=trim(OUTPUT_FILES)//'/output_write_profile.txt',status='unknown')

    write(IMAIN,*) '**************'
    write(IMAIN,*) 'xwrite_profile'
    write(IMAIN,*) '**************'

    ! command line when calling
    call get_command(arg)
    write(IMAIN,*) 'called by: ',trim(arg)
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! initializes
  ! default 2x2 degree grid
  delta_lat = INCR_lat
  delta_lon = INCR_lon

  ! start position
  initial_colat =  COLAT_0! 90.d0 - start_lat
  initial_lon = LON_0
  use_userdefined_lat = .false.

  ! loop range
  icolat_start = COLAT_istart
  icolat_end = COLAT_iend
  ilon_start = LON_istart
  ilon_end = LON_iend

  ! checks program arguments
  count = command_argument_count()
  if (count == 0 .or. (count /= 2 .and. count /= 4)) call usage()

  ! gets increments
  do i = 1,count
    call get_command_argument(i,arg)
    ! usage info
    if (len_trim(arg) == 0) then
      call usage()
    else
      select case(i)
      case (1)
        read(arg,*) delta_lat
        ! new loop range
        if (delta_lat == 0.d0) then
          icolat_end = 0
        else
          icolat_end = int((180.0 - initial_colat)/delta_lat)
        endif
      case (2)
        read(arg,*) delta_lon
        ! new loop range
        if (delta_lon == 0.d0) then
          ilon_end = 0
        else
          ilon_end = int((360.0 - initial_lon)/delta_lon)
        endif
      case (3)
        read(arg,*) lat
        ! initial colatitude
        initial_colat = 90.d0 - lat
        ! single profile only
        icolat_start = 0
        icolat_end = 0
        ! set flag
        use_userdefined_lat = .true.
      case (4)
        read(arg,*) lon
        ! initial longitude
        initial_lon = lon
        ! single profile only
        ilon_start = 0
        ilon_end = 0
      end select
    endif
  enddo
  call synchronize_all()

  ! user output
  if (myrank == 0) then
    print *,'write profile:'
    print *,'  initial   latitude/longitude = ',sngl(90.d0 - initial_colat),sngl(initial_lon), '(degrees)'
    print *,'          colatitude/longitude = ',sngl(initial_colat),sngl(initial_lon), '(degrees)'
    print *
    print *,'  colatitude loop start/end = ',icolat_start,icolat_end
    print *,'  increment lat = ',delta_lat
    print *
    print *,'  longitude  loop start/end = ',ilon_start,ilon_end
    print *,'  increment lon = ',delta_lon
    print *
  endif

  ! gets parameters and loads crustal/mantle model
  call write_profile_setup()

! < this part is from create_regions_mesh -> initialize_layers()
  ! see also: call define_all_layers_number_and_offset(NUMBER_OF_MESH_LAYERS,idummy)
  if (ONE_CRUST) then
    NUMBER_OF_MESH_LAYERS = MAX_NUMBER_OF_MESH_LAYERS - 1
  else
    NUMBER_OF_MESH_LAYERS = MAX_NUMBER_OF_MESH_LAYERS
  endif
!> end of part from create_regions_mesh -> initialize_layers()

  ! corrects user defined starting lat by ellipticity factor to have geocentric latitude
  ! we will loop over geocentric positions and convert back to geographic for getting topography in case
  if (use_userdefined_lat .and. ELLIPTICITY) then
    ! gets geocentric colat in rad
    call lat_2_geocentric_colat_dble(lat,theta,ELLIPTICITY)
    ! re-sets initial colatitude
    initial_colat = theta * RADIANS_TO_DEGREES
  endif

  ! loop over all theta (colatitude) and phis (longitude), every two degrees
  do i = icolat_start,icolat_end
    do j = ilon_start,ilon_end
      ! target location
      theta_degrees = initial_colat + i*delta_lat ! colatitude [0,180]
      phi_degrees   = initial_lon   + j*delta_lon ! longitude [0,360]

      ! checks limits
      if (theta_degrees < 0.d0) stop 'Error invalid colatitude < 0'
      if (theta_degrees > 180.d0) stop 'Error invalid colatitude > 180'
      if (phi_degrees < 0.d0) phi_degrees = phi_degrees + 360.d0
      if (phi_degrees > 360.d0) phi_degrees = phi_degrees - 360.d0

      ! loads corresponding GLL mesh
      if (MODEL_GLL) call load_GLL_mesh(theta_degrees,phi_degrees)

      ! file name
      write(outfile,'(a,i04.4,a,i04.4)') trim(OUTPUT_FILES)//'/CARDS_th',int(theta_degrees),'_ph',int(phi_degrees)
      if (myrank == 0) print *, 'file: ',trim(outfile)

      ! colat/lon in rad
      theta = theta_degrees * PI/180.0d0
      phi   = phi_degrees   * PI/180.0d0

      ! theta between 0 and PI, and phi between 0 and 2*PI
      call reduce(theta,phi)

      ! gets moho depth
      call write_profile_moho_depth(theta,phi,moho)

      ! gets topography elevation
      call write_profile_elevation(theta,phi,elevation)

      ! open output file
      open(unit=57,file=trim(outfile),status='unknown',iostat=ier)
      if (ier /= 0) stop 'Error opening file'

      ! header info
      if (MODEL_GLL) then
        write(57,'(a)') '# model         : '//trim(MODEL)//' ! instead plotting values of GLL_REFERENCE_MODEL (see constants.h)'
      else
        write(57,'(a)') '# model         : '//trim(MODEL)
      endif
      if (ELLIPTICITY) then
        ! output geocentric and geographic position
        ! converts geocentric theta/phi position to geographic lat/lon (in degrees)
        call thetaphi_2_geographic_latlon_dble(theta,phi,lat,lon,ELLIPTICITY)
        write(57,'(a,F10.4)') '# geographic point location: colatitude [degree] = ',sngl(90.0 - lat)
        write(57,'(a,F10.4)') '#            point location: longitude  [degree] = ',sngl(lon)
        write(57,'(a,F10.4)') '# geocentric point location: colatitude [degree] = ',sngl(theta_degrees)
        write(57,'(a,F10.4)') '#            point location: longitude  [degree] = ',sngl(phi_degrees)
      else
        ! geocentric and geographic lat/lon are the same
        write(57,'(a,F10.4)') '# point location: colatitude [degree] = ',sngl(theta_degrees)
        write(57,'(a,F10.4)') '# point location: longitude  [degree] = ',sngl(phi_degrees)
      endif
      if (CRUSTAL) then
        write(57,'(a,F10.4)') '# moho depth    : moho [km] = ',sngl(moho*(R_PLANET/1000.d0))
      else
        write(57,'(a,F10.4)') '# moho 1D depth : moho [km] =',sngl((1.d0-RMOHO/R_PLANET)*(R_PLANET/1000.d0))
      endif
      !write(57,'(a)') '#'
      if (TRANSVERSE_ISOTROPY) then
        write(57,'(a)') '# transverse isotropy'
      else
        write(57,'(a)') '# no transverse isotropy'
      endif
      if (TOPOGRAPHY) then
        write(57,'(a,F10.4)') '# surface topography: elevation [m] = ',sngl(elevation)
      else
        write(57,'(a)') '# no surface topography'
      endif
      if (OCEANS .or. ((.not. CRUSTAL) .and. (ROCEAN < R_PLANET))) then
        if (OCEANS) then
          write(57,'(a)') '# oceans'
        else
          ! 1D crustal models
          if (TOPOGRAPHY) then
            if (PLANET_TYPE == IPLANET_EARTH) then
              write(57,'(a)') '# oceans'
            else
              write(57,'(a)') '# no oceans, only topography'
            endif
          else
            write(57,'(a)') '# 1D reference oceans'
          endif
        endif
      else
        if (TOPOGRAPHY) then
          if (PLANET_TYPE == IPLANET_EARTH) then
            write(57,'(a)') '# oceans'
          else
            write(57,'(a)') '# no oceans, only topography'
          endif
        else
          write(57,'(a)') '# no oceans'
        endif
      endif
      if (ELLIPTICITY) then
        write(57,'(a)') '# ellipticity'
      else
        write(57,'(a)') '# no ellipticity'
      endif
      write(57,'(a)') '#'
      write(57,'(a)') '#radius(m) #rho #vpv #vsv #Qkappa #Qmu #vph #vsh #eta'

      ! keep track of line number to be able to write out locations of discontinuities at end
      iline = 0
      iline_icb = -1
      iline_cmb = -1
      iline_moho = -1
      iline_ocean = -1
      rmax_last = 0.0d0

      ! loop over all layers
      do ilayer = 1,NUMBER_OF_MESH_LAYERS

        ! loops from inner core to surface
        if (ilayer == 1) then
          rmin = 0.0d0
          rmax = rmins(NUMBER_OF_MESH_LAYERS-1)
          idoubling = IFLAG_INNER_CORE_NORMAL
        else
          rmin = rmins(NUMBER_OF_MESH_LAYERS-ilayer+1)
          rmax = rmaxs(NUMBER_OF_MESH_LAYERS-ilayer+1)
          idoubling = doubling_index(NUMBER_OF_MESH_LAYERS-ilayer+1)
        endif

        ! sets region code
        select case (idoubling)
        case (IFLAG_INNER_CORE_NORMAL,IFLAG_MIDDLE_CENTRAL_CUBE,IFLAG_BOTTOM_CENTRAL_CUBE, &
              IFLAG_TOP_CENTRAL_CUBE,IFLAG_IN_FICTITIOUS_CUBE)
          iregion_code = IREGION_INNER_CORE

        case (IFLAG_OUTER_CORE_NORMAL)
          iregion_code = IREGION_OUTER_CORE

        case (IFLAG_MANTLE_NORMAL,IFLAG_670_220,IFLAG_220_80,IFLAG_80_MOHO,IFLAG_CRUST)
          iregion_code = IREGION_CRUST_MANTLE

        case default
          print *,'Error idoubling ',idoubling
          stop 'Invalid idoubling value'
        end select

        !  make sure that the Moho discontinuity is at the real moho
        if (CRUSTAL) then
          ! checks rmin == RMOHO_FICTITIOUS_IN_MESHER/R_PLANET
          if (abs(rmin - RMOHO_FICTITIOUS_IN_MESHER/R_PLANET) < TOL_ZERO) rmin = 1.0d0 - moho
          ! checks rmax == RMOHO_FICTITIOUS_IN_MESHER/R_PLANET
          if (abs(rmax - RMOHO_FICTITIOUS_IN_MESHER/R_PLANET) < TOL_ZERO) rmax = 1.0d0 - moho
          !print *,'rmin == moho at line ',iline
        endif

        if (abs(rmin - rmax_last) < TOL_ZERO) then !!!! rmin == rmax_last: this means that we have just jumped between layers
          ! depth increment
          ! write values every 10 km in the deep earth and every 1 km in the shallow earth
          if (rmin > ((R_PLANET/1000.d0)-DELTA_HIRES_DEPTH)/(R_PLANET/1000.d0)) then
            delta = 1.0d0   ! in km
          else
            delta = 10.0d0  ! in km
          endif

          ! normalization
          delta = delta/(R_PLANET/1000.d0)

          ! sets maximum radius without ocean for 1D models
          if (((.not. CRUSTAL) .and. (ROCEAN < R_PLANET)) .and. (.not. TOPOGRAPHY)) then
            ! stops at ocean depth and adds last ocean layers explicitly
            if (abs(rmax - 1.0d0) < TOL_ZERO) rmax = ROCEAN/R_PLANET      ! rmax == 1.d0
          endif

          ! backup to detect jump between layers
          rmax_last = rmax

          ! number of iterations in increments of delta between rmin and rmax
          ! note: instead of (rmax - rmin), we add a factor (rmax * 0.999999 - rmin) to avoid getting an extra step
          !       in case the difference is an exact delta match, since we add +1 to nit to reach rmax
          nit = floor((rmax*0.9999999d0 - rmin)/delta) + 1

          ! debug
          !print *,'debug: write profile ilayer/iregion ',ilayer,iregion_code,'rmin/rmax',rmin,rmax,'delta',delta,'nit',nit

          do idep = 1,nit+1
            ! line counters
            ! inner core boundary
            if (abs(rmin - RICB/R_PLANET) < TOL_ZERO .and. idep == 1) iline_icb = iline  ! rmin == RICB/R_PLANET
            ! core mantle boundary
            if (abs(rmin - RCMB/R_PLANET) < TOL_ZERO .and. idep == 1) iline_cmb = iline  ! rmin == RCMB/R_PLANET
            ! moho
            if (CRUSTAL) then
              ! uses 3D crustal model (e.g. Crust2.0)
              if (abs(rmin - (1.0d0 - moho)) < TOL_ZERO .and. idep == 1) then  ! rmin == (1.0d0 - moho)
                iline_moho = iline
              endif
            else
              ! 1D crust from reference model
              if (abs(rmin - RMOHO/R_PLANET) < TOL_ZERO .and. idep == 1) iline_moho = iline     ! rmin == RMOHO/R_PLANET
            endif

            ! radius
            r = rmin + (idep-1)*delta

            ! debug
            !print *,'debug: write profile radius ',ilayer,iregion_code,'idep',idep,nit+1,'r',r,rmin,'delta',delta

            ! make sure we are within the right shell in PREM to honor discontinuities
            ! use small geometrical tolerance
            r_prem = r
            if (r < rmin*1.000001d0) r_prem = rmin*1.000001d0
            if (r > rmax*0.999999d0) r_prem = rmax*0.999999d0

            ! gets model properties (similar to get_model() routine)
            call write_profile_model_values(r,r_prem,theta,phi,iregion_code,idoubling,rmin,rmax, &
                                            rho,vpv,vph,vsv,vsh,eta_aniso,Qmu,Qkappa)

            ! make sure that the first and last point match the layer boundaries
            if (idep == 1) then
              r_prem = rmin
            else if (idep == nit+1) then
              r_prem = rmax
            endif

            ! adds topography stretching
            call write_profile_add_topography(r_prem,idoubling,elevation)

            ! ellipticity
            if (ELLIPTICITY) then
              ! adds ellipticity factor to radius
              call add_ellipticity_rtheta(r_prem,theta,nspl_ellip,rspl_ellip,ellipicity_spline,ellipicity_spline2)
            endif

            ! scale values read from routines back to true values
            scaleval = dsqrt(PI*GRAV*RHOAV)
            rho = rho*RHOAV/1000.0d0              ! [kg/m3]
            vpv  = vpv*R_PLANET*scaleval/1000.0d0  ! [m/s]
            vph  = vph*R_PLANET*scaleval/1000.0d0  ! [m/s]
            vsv  = vsv*R_PLANET*scaleval/1000.0d0  ! [m/s]
            vsh  = vsh*R_PLANET*scaleval/1000.0d0  ! [m/s]

            ! finally write the values obtained at the given depth to file
            str_info = ''
            if (iline == iline_icb) str_info = ' # ICB'
            if (iline == iline_cmb) str_info = ' # CMB'
            if (iline == iline_moho) str_info = ' # moho'

            write(57,'(F10.0,7F12.2,F12.5,a)') &
              sngl(r_prem*R_PLANET),sngl(rho*1000.d0),sngl(vpv*1000.d0),sngl(vsv*1000.d0), &
              sngl(Qkappa),sngl(Qmu),sngl(vph*1000.d0),sngl(vsh*1000.d0),sngl(eta_aniso),trim(str_info)

            ! line counter
            iline = iline + 1

            ! debug
            !write(*,'(i3,11F12.4)') &
            ! iline,sngl(rmin*(R_PLANET/1000.d0)),sngl(rmax*(R_PLANET/1000.d0)), &
            ! sngl(r_prem*(R_PLANET/1000.d0)),sngl(r*(R_PLANET/1000.d0)), &
            ! sngl(vpv),sngl(vph),sngl(vsv),sngl(vsh),sngl(rho),sngl(eta_aniso),sngl(Qmu)

          enddo !idep
        endif !rmin == rmax_last
      enddo !ilayer

      ! adds ocean layer profile
      call write_profile_ocean(r_prem,theta,elevation,iline,iline_ocean)

      ! discontinuities info
      write(57,'(a)') '# discontinuities info'
      write(57,'(a)') '# line number (not counting comment lines #..) for'
      write(57,'(a,i5)') '# last line = ',iline
      write(57,'(a,i5)') '# ICB       = ',iline_icb
      write(57,'(a,i5)') '# CMB       = ',iline_cmb
      write(57,'(a,i5)') '# Moho      = ',iline_moho
      if (iline_ocean > 0) write(57,'(a,i5)') '# Ocean     = ',iline_ocean

    enddo !sum over phi
  enddo !sum over theta

  ! closes main output file
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'see profiles CARDS_** written in directory: ',trim(OUTPUT_FILES)
    write(IMAIN,*)
    write(IMAIN,*) 'End of xwrite_profile'
    write(IMAIN,*)
    call flush_IMAIN()
    close(IMAIN)

    print *
    if (IMAIN /= ISTANDARD_OUTPUT) &
      print *,'output written to: ',trim(OUTPUT_FILES)//'/output_write_profile.txt'
  endif

  call finalize_mpi()

  end program xwrite_profile

!
!-------------------------------------------------------------------------------------
!

  subroutine usage()

  use constants, only: myrank

  implicit none

  if (myrank == 0) then
    print *,'Usage: xwrite_profile [dlat] [dlon] [lat0] [lon0]'
    print *,'   with'
    print *,'      dlat dlon         - (optional) increments for latitude/longitude (in degrees, by default 2x2 degrees)'
    print *,'      lat0 lon0         - (optional) latitude/longitude (in degrees) for single profile output'
    print *
  endif
  stop ' Reenter command line options'

  end subroutine usage

!
!-------------------------------------------------------------------------------------
!

  subroutine write_profile_setup()

  use constants, only: IMAIN,myrank,N_SLS,NGLLX,NGLLY,NGLLZ,SUPPRESS_MOHO_STRETCHING
  use meshfem_models_par
  use shared_parameters

  implicit none

! < this part is from meshfem3D and setup_model (read parfile and output info)
  if (myrank == 0) print *,'reading parameter file...'

  ! read the parameter file and compute additional parameters
  call read_compute_parameters()

  ! user output
  if (myrank == 0) then
    ! model
    print *
    print *, 'model: ',trim(MODEL)
    print *

    ! file output
    write(IMAIN,*) 'model: ',trim(MODEL)
    ! mesh streching
    if (CRUSTAL .and. CASE_3D) then
      if (REGIONAL_MOHO_MESH) then
        if (SUPPRESS_MOHO_STRETCHING) then
          write(IMAIN,*) '  no element stretching for 3-D moho surface'
        else
          write(IMAIN,*) '  incorporating element stretching for regional 3-D moho surface'
        endif
      else
        if (SUPPRESS_MOHO_STRETCHING .or. (.not. TOPOGRAPHY)) then
          write(IMAIN,*) '  no element stretching for 3-D moho surface'
        else
          write(IMAIN,*) '  incorporating element stretching for 3-D moho surface'
        endif
      endif
    else
      write(IMAIN,*) '  no element stretching for 3-D moho surface'
    endif
    write(IMAIN,*)
    ! model user parameters
    if (OCEANS) then
      write(IMAIN,*) '  incorporating the oceans using equivalent load'
    else
      write(IMAIN,*) '  no oceans'
    endif
    if (ELLIPTICITY) then
      write(IMAIN,*) '  incorporating ellipticity'
    else
      write(IMAIN,*) '  no ellipticity'
    endif
    if (TOPOGRAPHY) then
      write(IMAIN,*) '  incorporating surface topography'
    else
      write(IMAIN,*) '  no surface topography'
    endif
    if (GRAVITY) then
      write(IMAIN,*) '  incorporating self-gravitation (Cowling approximation)'
    else
      write(IMAIN,*) '  no self-gravitation'
    endif
    if (ROTATION) then
      write(IMAIN,*) '  incorporating rotation'
    else
      write(IMAIN,*) '  no rotation'
    endif
    if (ATTENUATION) then
      write(IMAIN,*) '  incorporating attenuation using ',N_SLS,' standard linear solids'
      if (ATTENUATION_3D) write(IMAIN,*) '  using 3D attenuation model'
    else
      write(IMAIN,*) '  no attenuation'
    endif
    write(IMAIN,*)
    ! model mesh parameters
    if (MODEL_3D_MANTLE_PERTUBATIONS) then
      write(IMAIN,*) '  incorporating 3-D lateral variations in the mantle'
    else
      write(IMAIN,*) '  no 3-D lateral variations in the mantle'
    endif
    if (HETEROGEN_3D_MANTLE) then
      write(IMAIN,*) '  incorporating heterogeneities in the mantle'
    else
      write(IMAIN,*) '  no heterogeneities in the mantle'
    endif
    if (CRUSTAL) then
      write(IMAIN,*) '  incorporating crustal variations'
    else
      write(IMAIN,*) '  no crustal variations'
    endif
    if (ONE_CRUST) then
      write(IMAIN,*) '  using one layer only in crust'
    else
      write(IMAIN,*) '  using unmodified 1D crustal model with two layers'
    endif
    if (TRANSVERSE_ISOTROPY) then
      write(IMAIN,*) '  incorporating transverse isotropy'
    else
      write(IMAIN,*) '  no transverse isotropy'
    endif
    if (ANISOTROPIC_INNER_CORE) then
      write(IMAIN,*) '  incorporating anisotropic inner core'
    else
      write(IMAIN,*) '  no inner-core anisotropy'
    endif
    if (ANISOTROPIC_3D_MANTLE) then
      write(IMAIN,*) '  incorporating anisotropic mantle'
    else
      write(IMAIN,*) '  no general mantle anisotropy'
    endif

    ! model mesh setup
    write(IMAIN,*)
    write(IMAIN,*) 'mesh setup: '
    write(IMAIN,*) '  NCHUNKS            : ',NCHUNKS
    write(IMAIN,*) '  NSPEC crust/mantle : ',NSPEC_REGIONS(IREGION_CRUST_MANTLE)
    write(IMAIN,*) '  NSPEC outer core   : ',NSPEC_REGIONS(IREGION_OUTER_CORE)
    write(IMAIN,*) '  NSPEC inner core   : ',NSPEC_REGIONS(IREGION_INNER_CORE)
    write(IMAIN,*) '  NGLLX/NGLLY/NGLLZ  : ',NGLLX,NGLLY,NGLLZ
    write(IMAIN,*)
    write(IMAIN,*)

    call flush_IMAIN()
  endif


  if (myrank == 0) print *,'reading model...'

  ! loads 3D models
  call meshfem3D_models_broadcast()

  ! user output
  if (myrank == 0) print *

!> end part from meshfem3D

  end subroutine write_profile_setup

!
!-------------------------------------------------------------------------------------
!


  subroutine write_profile_moho_depth(theta,phi,moho)

  use constants, only: myrank
  use meshfem_models_par
  use shared_parameters

  implicit none

  double precision,intent(in) :: theta,phi
  double precision,intent(inout) :: moho

  ! local parameters
  double precision :: r,sediment
  double precision :: rho,vpv,vph,vsv,vsh,eta_aniso
  double precision :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  ! moho depth
  ! read crustal models and topo models as they are needed to modify the depths of the discontinuities
  if (CRUSTAL) then
    ! fixes depth to r=-7km
    r = 1.0d0 - 7.0d0/(R_PLANET/1000.d0)

    !if (myrank == 0) print *, 'debug: r,theta,phi = ',r,90.0-theta*180./PI,phi*180./PI

    ! gets moho depth from crustal model
    call meshfem3D_models_get3Dcrust_val(IREGION_CRUST_MANTLE,r,theta,phi, &
                                         vpv,vph,vsv,vsh,rho,eta_aniso, &
                                         c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                         c33,c34,c35,c36,c44,c45,c46,c55,c56,c66, &
                                         .true.,moho,sediment)

    !if (myrank == 0) print *, '  crustal values:',vpv,vph,vsv,vsh,rho,eta_aniso
    if (myrank == 0) print *, '  moho depth [km]    :',sngl(moho*(R_PLANET/1000.d0))
                              !,'radius (in km):',sngl((1.0d0-moho)*(R_PLANET/1000.d0))
    if (myrank == 0) print *, '  sediment depth [km]:',sngl(sediment*(R_PLANET/1000.d0))
  else
    moho = 0.d0
    ! 1D crust from reference model
    if (myrank == 0) print *, '  moho 1D depth [km]:',sngl((1.d0-RMOHO/R_PLANET)*(R_PLANET/1000.d0))
  endif

  end subroutine write_profile_moho_depth

!
!-------------------------------------------------------------------------------------
!

  subroutine write_profile_elevation(theta,phi,elevation)

  use constants, only: THICKNESS_OCEANS_PREM,myrank
  use meshfem_models_par
  use shared_parameters

  implicit none

  double precision,intent(in) :: theta,phi
  double precision,intent(out) :: elevation

  ! local parameters
  double precision :: lat,lon

  ! initializes
  elevation = 0.d0

  ! topography elevation
  if (TOPOGRAPHY .or. OCEANS) then
    if (TOPOGRAPHY) then
      ! converts geocentric theta/phi position to geographic lat/lon (in degrees)
      call thetaphi_2_geographic_latlon_dble(theta,phi,lat,lon,ELLIPTICITY)

      ! longitude range [-180,180]
      if (lon > 180.0d0) lon = lon-360.0d0

      !if (myrank == 0) print *,'  get_topo_bathy(lat,lon,elevation,ibathy_topo',lat,lon,elevation

      call get_topo_bathy(lat,lon,elevation,ibathy_topo)

      !if (myrank == 0) print *,'  get_topo_bathy(lat,lon,elevation,ibathy_topo',lat,lon,elevation
      if (myrank == 0) print *,'  elevation [m]:',sngl(elevation) !, &
                               !'surface radius (in km):',sngl((1.0d0 + elevation/R_PLANET)*(R_PLANET/1000.d0))

    else
      ! oceans selected only, without topography
      if ((.not. CRUSTAL) .and. (ROCEAN < R_PLANET)) then
        ! if 1D Earth, use oceans of constant thickness everywhere (see create_mass_matrices_ocean_load())
        elevation = - THICKNESS_OCEANS_PREM * R_PLANET ! in m
        if (myrank == 0) print *,'  ocean PREM [m]:',sngl(elevation)
      else
        ! 3D crustal model without topography, not using an ocean
        elevation = 0.d0
        if (myrank == 0) print *,'  elevation [m]:',sngl(elevation)
      endif
    endif
  else
    elevation = 0.d0
    if (myrank == 0) print *,'  elevation [m]:',sngl(elevation)
  endif

  end subroutine write_profile_elevation

!
!-------------------------------------------------------------------------------------
!

  subroutine write_profile_model_values(r,r_prem,theta,phi,iregion_code,idoubling,rmin,rmax, &
                                        rho,vpv,vph,vsv,vsh,eta_aniso,Qmu,Qkappa)

  use constants, only: myrank,N_SLS,TINYVAL
  use meshfem_models_par
  use shared_parameters

  implicit none

  double precision,intent(in) :: r,r_prem,theta,phi
  integer,intent(in) :: iregion_code,idoubling
  double precision,intent(in) :: rmin,rmax
  double precision,intent(out) :: rho,vpv,vph,vsv,vsh,eta_aniso,Qmu,Qkappa

  ! local parameters
  double precision :: moho,sediment
  double precision :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66
  ! Attenuation values
  double precision, dimension(N_SLS) :: tau_s, tau_e
  ! Parameter used to decide whether this element is in the crust or not
  logical :: elem_in_crust,elem_in_mantle
  integer :: i,j,k,ispec

  ! initializes values
  rho = 0.d0
  vpv = 0.d0
  vph = 0.d0
  vsv = 0.d0
  vsh = 0.d0
  eta_aniso = 1.d0
  c11 = 0.d0
  c12 = 0.d0
  c13 = 0.d0
  c14 = 0.d0
  c15 = 0.d0
  c16 = 0.d0
  c22 = 0.d0
  c23 = 0.d0
  c24 = 0.d0
  c25 = 0.d0
  c26 = 0.d0
  c33 = 0.d0
  c34 = 0.d0
  c35 = 0.d0
  c36 = 0.d0
  c44 = 0.d0
  c45 = 0.d0
  c46 = 0.d0
  c55 = 0.d0
  c56 = 0.d0
  c66 = 0.d0
  Qmu = 0.d0
  Qkappa = 0.d0 ! not used, not stored so far...
  tau_e(:) = 0.d0
  tau_s(:) = 0.d0
  moho = 0.d0
  sediment = 0.d0

  i = 1; j = 1; k = 1; ispec = 1 ! dummy GLL point

  ! do not force the element/point to be in the crust
  elem_in_crust = .false.

! < start GET_MODEL
  ! checks r_prem,rmin/rmax and assigned idoubling
  call get_model_check_idoubling(r_prem,theta,phi,rmin,rmax,idoubling, &
                                 RICB,RCMB,RTOPDDOUBLEPRIME, &
                                 R220,R670)

  ! gets reference model values: rho,vpv,vph,vsv,vsh and eta_aniso
  call meshfem3D_models_get1D_val(iregion_code,idoubling, &
                                  r_prem,rho,vpv,vph,vsv,vsh,eta_aniso, &
                                  Qkappa,Qmu,RICB,RCMB, &
                                  RTOPDDOUBLEPRIME,R80,R120,R220,R400,R670,R771, &
                                  RMOHO,RMIDDLE_CRUST)

  ! gets the 3-D model parameters for the mantle
  call meshfem3D_models_get3Dmntl_val(iregion_code,r_prem,rho, &
                                      vpv,vph,vsv,vsh,eta_aniso, &
                                      RCMB,RMOHO, &
                                      r,theta,phi, &
                                      c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                      c33,c34,c35,c36,c44,c45,c46,c55,c56,c66, &
                                      ispec,i,j,k)

  ! gets the 3-D crustal model
  if (CRUSTAL) then
    ! taking r_prem to be within layer (slight difference to get_model() routine)
    ! avoid calling the 3Dcrust routine for points below moho
    if (r_prem < R_DEEPEST_CRUST) then
      elem_in_mantle = .true.
    else
      elem_in_mantle = .false.
    endif
    if (.not. elem_in_mantle) then
      call meshfem3D_models_get3Dcrust_val(iregion_code,r_prem,theta,phi, &
                                           vpv,vph,vsv,vsh,rho,eta_aniso, &
                                           c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                           c33,c34,c35,c36,c44,c45,c46,c55,c56,c66, &
                                           elem_in_crust,moho,sediment)

      !print *,'crustal: ',r_prem,R_DEEPEST_CRUST,elem_in_crust,vpv,vph,vsv,vsh,eta_aniso
    endif
  endif

  !!! VH  commented out following two lines from get_model
  !! overwrites with tomographic model values (from iteration step) here, given at all GLL points
  !call meshfem3D_models_impose_val(iregion_code,r,theta,phi,ispec,i,j,k, &
  !                                 vpv,vph,vsv,vsh,rho,eta_aniso, &
  !                                 c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
  !                                 c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)
  !
  ! note: for GLL models, we would need to have the GLL point index (ispec,i,j,k) instead of location (r_prem,theta,phi).
  !       the mesh indexing however is only available after running the full mesher.
  !
  !       todo: a workaround would be to run first mesher, then read in the database files and search for the closest GLL
  !             point to a given location (x,y,z).
  !
  !       for now, we will output the values from the initial reference model used by the GLL model,
  !       for example s362ani instead of GLL as choosen by GLL_REFERENCE_MODEL in constants.h

  ! checks vpv: if close to zero then there is probably an error
  if (vpv < TINYVAL) then
    print *,'error vpv: ',vpv,vph,vsv,vsh,rho
    print *,'radius:',r*(R_PLANET/1000.d0)
    call exit_mpi(myrank,'Error get model values in write_profile')
  endif

  !> Hejun
  ! New Attenuation assignment
  ! Define 3D and 1D Attenuation after moho stretch
  ! and before TOPOGRAPHY/ELLIPCITY
  !
  !note:  only Qmu attenuation considered, Qkappa attenuation not used so far...
  if (ATTENUATION) then
    call meshfem3D_models_getatten_val(idoubling,r_prem,theta,phi, &
                                       ispec, i, j, k, &
                                       tau_e,tau_s, &
                                       moho,Qmu,Qkappa,elem_in_crust)
  endif
!> end GET_MODEL

  end subroutine write_profile_model_values

!
!-------------------------------------------------------------------------------------
!


  subroutine write_profile_add_topography(r_prem,idoubling,elevation)

  use constants, only: IFLAG_CRUST,IFLAG_80_MOHO,IFLAG_220_80,R_UNIT_SPHERE
  use shared_parameters
  use meshfem_models_par

  implicit none

  double precision,intent(inout) :: r_prem
  double precision,intent(in) :: elevation
  integer,intent(in) :: idoubling

  ! local parameters
  ! for stretching
  double precision :: gamma

! < start add_topography
  if (TOPOGRAPHY) then
    if (idoubling == IFLAG_CRUST .or. &
        idoubling == IFLAG_80_MOHO .or. &
        idoubling == IFLAG_220_80) then
      ! adds mesh stretching
      !print *, 'adding topography.  elevation: ',elevation
      gamma = (r_prem - R220/R_PLANET) / (R_UNIT_SPHERE - R220/R_PLANET)
      if (gamma < -0.02 .or. gamma > 1.02) print *, 'incorrect value of gamma for topograpy'

      ! stretching due to elevation
      !print *,'rprem before: ',r_prem*R_PLANET
      r_prem = r_prem*(ONE + gamma * (elevation/R_PLANET) /r_prem)
      !print *,'r_prem after: ',r_prem*R_PLANET
    endif
  endif
!> end add_topography

  end subroutine write_profile_add_topography

!
!-------------------------------------------------------------------------------------
!


  subroutine write_profile_ocean(r_prem,theta,elevation,iline,iline_ocean)

  use constants, only: MINIMUM_THICKNESS_3D_OCEANS
  use meshfem_models_par
  use shared_parameters

  implicit none

  double precision,intent(in) :: r_prem,theta,elevation
  integer,intent(inout) :: iline,iline_ocean

  ! local parameters
  ! depth increment
  double precision :: delta

  ! ocean
  integer :: ilayers_ocean,nlayers_ocean
  double precision :: r_ocean

! < ocean
  ! This part adds the ocean to profile where needed
  if ( (TOPOGRAPHY .or. OCEANS) .or. ((.not. CRUSTAL) .and. (ROCEAN < R_PLANET)) ) then

    ! note: For 1D models (with 1D reference crust), CRUSTAL is set to .false.
    !
    !       The meshing by meshfem3D is not including ocean values, instead it will use upper crust velocities.
    !       This is mainly due to difficulties to mesh elements in shallow oceans. To account for the effect of oceans
    !       on waveforms, the ocean loading is approximated using the (actual) ocean thickness.
    !
    !       Here for the models, we want to output a profile which includes the (original) ocean layer.
    !       When TOPOGRAPHY is selected, the actual ocean thickness will be considered as in the meshing part
    !       with a minimum thickness threshold.
    !
    !       For cases where OCEANS is selected to include the ocean loading effect on waveforms, a corresponding mass matrix
    !       will be created which includes the actual ocean thickness if TOPOGRAPHY is selected or
    !       a default PREM ocean thickness (3km).

    ! initializes ocean profile layers
    nlayers_ocean = 0
    r_ocean = 1.d0

    ! depth increment
    delta = 500.d0 ! in m

    ! determines final ocean layers (nlayers_ocean)
    if (OCEANS .and. elevation < - MINIMUM_THICKNESS_3D_OCEANS) then
      if (TOPOGRAPHY) then
        ! uses actual bathymetry value and fills up ocean to surface
        nlayers_ocean = floor(-elevation/delta)
      else
        ! no topography on mesh
        if ((.not. CRUSTAL) .and. (ROCEAN < R_PLANET)) then
          ! 1D model without topography, uses default ocean layer (R_PLANET - ROCEAN)
          nlayers_ocean = floor((R_PLANET - ROCEAN)/delta)
        else
          ! 3D crustal model, without topography
          ! no oceans
          nlayers_ocean = -1
        endif
      endif
    else
      ! case for 1D (reference crust) models which include an ocean: ((.not. CRUSTAL) .and. (ROCEAN < R_PLANET)) then
      if (TOPOGRAPHY) then
        if (elevation < - MINIMUM_THICKNESS_3D_OCEANS) then
          ! uses actual bathymetry value and fills up ocean to surface
          nlayers_ocean = floor(-elevation/delta) ! this means ocean depth must be bigger than 500m (delta) for output
        else
          ! elevation above minimum ocean depth
          nlayers_ocean = -1
        endif
      else
        if ((.not. CRUSTAL) .and. (ROCEAN < R_PLANET)) then
          ! 1D models without topography, takes reference ocean depth
          nlayers_ocean = floor((R_PLANET - ROCEAN)/delta)
        else
          ! 3D crustal models without topography will have no ocean
          nlayers_ocean = -1
        endif
      endif
    endif

    ! Mars and Moon models have no ocean layers
    if (PLANET_TYPE /= IPLANET_EARTH) nlayers_ocean = 0

    ! adds ocean layer
    if (nlayers_ocean > 0) then
      ! ocean line
      iline_ocean = iline

      ! points within ocean
      do ilayers_ocean = 0,nlayers_ocean
        ! radius
        r_ocean = r_prem + ilayers_ocean * delta/R_PLANET

        ! ellipticity
        if (ELLIPTICITY) then
          ! adds ellipticity factor to radius
          call add_ellipticity_rtheta(r_ocean,theta,nspl_ellip,rspl_ellip,ellipicity_spline,ellipicity_spline2)
        endif

        ! ocean properties (salt water parameters from PREM)
        if (ilayers_ocean == 0) then
          ! line with section info
          write(57,'(F10.0,7F12.2,F12.5,a)') sngl(r_ocean*R_PLANET),1020.0,1450.,0.0,57822.5,0.0,1450.0,0.0,1.0,' # ocean'
        else
          write(57,'(F10.0,7F12.2,F12.5)') sngl(r_ocean*R_PLANET),1020.0,1450.,0.0,57822.5,0.0,1450.0,0.0,1.0
        endif
        ! line counter
        iline = iline + 1
      enddo
      ! at surface
      if (r_ocean < 1.d0) then
        ! last line exactly at earth surface
        r_ocean = 1.d0
        ! ellipticity
        if (ELLIPTICITY) then
          ! adds ellipticity factor to radius
          call add_ellipticity_rtheta(r_ocean,theta,nspl_ellip,rspl_ellip,ellipicity_spline,ellipicity_spline2)
        endif
        ! last line
        write(57,'(F10.0,7F12.2,F12.5)') sngl(r_ocean*R_PLANET),1020.0,1450.,0.0,57822.5,0.0,1450.,0.0,1.0
        ! line counter
        iline = iline + 1
      endif
    endif
  endif
!> end adding ocean

  end subroutine write_profile_ocean

!
!-------------------------------------------------------------------------------------
!


  subroutine load_GLL_mesh(theta_degrees,phi_degrees)

! loads the GLL mesh slice which contains the colat/lon position

  use constants, only: CUSTOM_REAL,myrank

  use shared_parameters, only: NCHUNKS,NPROC_XI,NPROCTOT,ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES, &
    ADIOS_FOR_MODELS,TRANSVERSE_ISOTROPY

  use model_gll_par

  implicit none

  double precision,intent(in) :: theta_degrees,phi_degrees

  ! local parameters
  real(kind=CUSTOM_REAL) :: lat,lon,xi,eta,xi_width,eta_width
  integer :: chunk
  integer :: slice_proc
  integer,save :: slice_proc_old = 0

  ! converts to lat/lon
  lat = real(90.d0 - theta_degrees,kind=CUSTOM_REAL)
  lon = real(phi_degrees,kind=CUSTOM_REAL)
  ! limits value range
  if (lat < -90.0) lat = -90.0
  if (lat > 90.0) lat = 90.0
  if (lon < -180.0) lon = lon + 360.0
  if (lon > 180.0) lon = lon - 360.0

  call get_latlon_chunk_location(lat,lon,NCHUNKS,chunk,xi,eta)

  ! chunk width
  if (NCHUNKS == 6) then
    xi_width = 90.0
    eta_width = 90.0
  else
    xi_width = real(ANGULAR_WIDTH_XI_IN_DEGREES,kind=CUSTOM_REAL)
    eta_width = real(ANGULAR_WIDTH_ETA_IN_DEGREES,kind=CUSTOM_REAL)
  endif

  ! converts chunk width to radians
  call get_process_slice_number(chunk,xi,eta,NPROC_XI,xi_width,eta_width,slice_proc)

  ! checks
  if (slice_proc < 0 .or. slice_proc > NPROCTOT) then
    print *,"Error invalid slice number found ",slice_proc," with total number of slices ",NPROCTOT
    stop 'Error invalid slice_proc'
  endif

  ! user output
  if (myrank == 0) then
    print *
    print *,'GLL model: location lat/lon = ',lat,'/',lon,' has mesh slice number = ',slice_proc
  endif

  ! reads in new model files in case
  if (slice_proc /= slice_proc_old) then
    ! user output
    if (myrank == 0) then
      print *,'GLL model: reloading mesh slice ',slice_proc
    endif
    ! loads GLL model
    if (ADIOS_FOR_MODELS) then
      call read_gll_model_adios(slice_proc)
    else
      call read_gll_model(slice_proc)
    endif
    ! saves slice number
    slice_proc_old = slice_proc
  endif

  ! min/max values output
  if (myrank == 0) then
    if (.not. TRANSVERSE_ISOTROPY) then
      ! isotropic model
      print *,'  vs new min/max: ',minval(MGLL_V%vs_new),maxval(MGLL_V%vs_new)
      print *,'  vp new min/max: ',minval(MGLL_V%vp_new),maxval(MGLL_V%vp_new)
    else
      ! transverse isotropic model
      print *,'  vsv new min/max: ',minval(MGLL_V%vsv_new),maxval(MGLL_V%vsv_new)
      print *,'  vsh new min/max: ',minval(MGLL_V%vsh_new),maxval(MGLL_V%vsh_new)
      print *,'  vpv new min/max: ',minval(MGLL_V%vpv_new),maxval(MGLL_V%vpv_new)
      print *,'  vph new min/max: ',minval(MGLL_V%vph_new),maxval(MGLL_V%vph_new)
      print *,'  eta new min/max: ',minval(MGLL_V%eta_new),maxval(MGLL_V%eta_new)
    endif
    print *
  endif

  end subroutine load_GLL_mesh


!
!-------------------------------------------------------------------------------------
!

  subroutine get_latlon_chunk_location(lat,lon,nchunks,chunk,xi,eta)

! returns the chunk number and process number xi/eta for a given location lat/lon (degrees)

  use constants, only: NDIM,PI,CUSTOM_REAL
  use shared_parameters, only: CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH

  implicit none

  real(kind=CUSTOM_REAL),intent(in) :: lat,lon
  integer,intent(in) :: nchunks
  integer,intent(out) :: chunk
  real(kind=CUSTOM_REAL),intent(out) :: xi,eta

  ! local parameters
  real(kind=CUSTOM_REAL) :: r,theta,phi,x,y,z,xik,etak
  real(kind=CUSTOM_REAL) :: xn,yn,zn
  double precision, dimension(NDIM,NDIM) :: rotation_matrix
  integer :: k

  ! converts lat/lon from degrees to radians (colatitute/longitude)
  theta = real((90.d0 - lat)/180.d0 * PI,kind=CUSTOM_REAL)
  phi = real(lon/180.d0 * PI,kind=CUSTOM_REAL)

  ! converts (r,theta,phi) to (x,y,z) on unit sphere
  r = 1.0_CUSTOM_REAL
  call rthetaphi_2_xyz(x,y,z,r,theta,phi)

  if (nchunks /= 6) then
    ! compute rotation matrix from Euler angles
    call euler_angles(rotation_matrix,CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH)

    xn = real(x * rotation_matrix(1,1) + y * rotation_matrix(2,1) + z * rotation_matrix(3,1),kind=CUSTOM_REAL)
    yn = real(x * rotation_matrix(1,2) + y * rotation_matrix(2,2) + z * rotation_matrix(3,2),kind=CUSTOM_REAL)
    zn = real(x * rotation_matrix(1,3) + y * rotation_matrix(2,3) + z * rotation_matrix(3,3),kind=CUSTOM_REAL)

    x = xn; y = yn; z = zn
  endif

  ! loops over chunks
  chunk = 0
  xi = 10.d0
  eta = 10.d0

  do k = 1,nchunks
    call chunk_map(k,x,y,z,xik,etak)

    ! checks location
    if (abs(xik) <= PI/4.0 .and. abs(etak) <= PI/4.0) then
      chunk = k
      xi = xik
      eta = etak
      !print *,"chunk ",k," contains slice with point"
      !print *,"   xi = ",xi," eta = ",eta
    endif

    ! mininum
    if (abs(xik) < xi .and. abs(etak) < eta) then
      chunk = k
      xi = xik
      eta = etak
    endif
  enddo

  if (nchunks /= 6 .and. (chunk > nchunks .or. abs(xi) > PI/2.0 .or. abs(eta) > PI/2.0)) then
    print *,'Check if point lat/lon = ', lat, '/', lon,' is within the region or not'
    stop 'Error invalid location for 1-chunk mesh'
  endif

  ! checks if we found a chunk
  if (chunk == 0) then
    print *,"no chunk found for location ",lat,"/",lon," , please check input..."
    stop 'Error chunk not found, invalid lat/lon'
  endif

  end subroutine get_latlon_chunk_location

!
!-------------------------------------------------------------------------------------
!
! not used any further... left here for reference
!
!  subroutine tp2xyz(th,ph,x,y,z)
!
!  use constants, only: CUSTOM_REAL
!
!  implicit none
!
!  real(CUSTOM_REAL),intent(in) :: th,ph
!  real(CUSTOM_REAL),intent(out) :: x,y,z
!
!  ! convert (th,ph) to (x,y,z) on unit sphere
!  x = sin(th) * cos(ph)
!  y = sin(th) * sin(ph)
!  z = cos(th)
!
!  end subroutine tp2xyz
!
!
!-------------------------------------------------------------------------------------
!
! not used any further... left here as a reference
!
!  subroutine rotmat(clon,clat,grot,rotation_matrix)
!
!  ! this function calculate the 3x3 rotation matrix from the AB chunk
!  ! frame to the actual frame defined by (clon,clat,grot)
!
!  use constants, only: CUSTOM_REAL,PI
!
!  implicit none
!
!  real(kind=CUSTOM_REAL),intent(in) :: clon,clat,grot
!  real(kind=CUSTOM_REAL),intent(out) :: rotation_matrix(3,3)
!
!  ! local parameters
!  real(kind=CUSTOM_REAL) :: alpha, beta, gamma, sina, cosa, sinb, cosb, sing, cosg
!
!  ! compute colatitude and longitude
!  alpha = clon
!  beta = PI/2.0 - clat
!  gamma = grot
!
!  sina = sin(alpha)
!  cosa = cos(alpha)
!  sinb = sin(beta)
!  cosb = cos(beta)
!  sing = sin(gamma)
!  cosg = cos(gamma)
!
!  ! define rotation matrix
!  rotation_matrix(1,1) = cosg*cosb*cosa-sing*sina
!  rotation_matrix(1,2) = -sing*cosb*cosa-cosg*sina
!  rotation_matrix(1,3) = sinb*cosa
!  rotation_matrix(2,1) = cosg*cosb*sina+sing*cosa
!  rotation_matrix(2,2) = -sing*cosb*sina+cosg*cosa
!  rotation_matrix(2,3) = sinb*sina
!  rotation_matrix(3,1) = -cosg*sinb
!  rotation_matrix(3,2) = sing*sinb
!  rotation_matrix(3,3) = cosb
!
!  end subroutine rotmat
!
!
!-------------------------------------------------------------------------------------
!

  subroutine chunk_map(k,xx,yy,zz,xi,eta)

  ! this routine gets the xi,eta for (xx,yy,zz)
  ! point under the k'th chunk coordinate
  ! transformation

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: k
  real(kind=CUSTOM_REAL),intent(in) :: xx, yy, zz
  real(kind=CUSTOM_REAL),intent(out) :: xi, eta

  ! local parameters
  real(kind=CUSTOM_REAL) :: x, y, z
  real(kind=CUSTOM_REAL), parameter :: EPS = 1e-6

  x = xx; y = yy; z = zz

  if (x < EPS .and. x >= 0) x = EPS
  if (x > -EPS .and. x < 0) x = -EPS
  if (y < EPS .and. y >= 0) y = EPS
  if (y > -EPS .and. y < 0) y = -EPS
  if (z < EPS .and. z >= 0) z = EPS
  if (z > -EPS .and. z < 0) z = -EPS

  if (k == 1) then
    ! CHUNK_AB
    xi = atan(y/z); eta = atan(-x/z)
    if (z < 0)  xi = 10.0
  else if (k == 2) then
    ! CHUNK_AC
    xi = atan(-z/y); eta = atan(x/y)
    if (y > 0) xi = 10.0
  else if (k == 3) then
    ! CHUNK_BC
    xi = atan(-z/x); eta = atan(-y/x)
    if (x > 0) xi = 10.0
  else if (k == 4) then
    ! CHUNK_AC'
    xi = atan(-z/y); eta = atan(-x/y)
    if (y < 0) xi = 10.0
  else if (k == 5) then
    ! CHUNK_BC'
    xi = atan(z/x); eta = atan(-y/x)
    if (x < 0) xi = 10.0
  else if (k == 6) then
    ! CHUNK_AB'
    xi = atan(y/z); eta = atan(x/z)
    if (z > 0)  xi = 10.0
  else
    stop 'chunk number k < 6'
  endif

  end subroutine chunk_map

!
!-------------------------------------------------------------------------------------
!

  subroutine get_process_slice_number(chunk,xi,eta,nproc_xi,xi_width,eta_width,slice_proc)

! returns slice number given chunk and chunk-location xi/eta,
! using a chunk width xi_width/eta_width and a number of nprocesses per chunk nproc_xi

  use constants, only: PI,CUSTOM_REAL

  implicit none

  integer,intent(in) :: chunk,nproc_xi
  real(kind=CUSTOM_REAL),intent(in) :: xi,eta
  real(kind=CUSTOM_REAL),intent(in) :: xi_width,eta_width

  integer,intent(out) :: slice_proc

  ! local parameters
  real(kind=CUSTOM_REAL) :: xi1,eta1,xi_width_rad,eta_width_rad
  integer :: proc_xi,proc_eta

  ! converts chunk width to radians
  xi_width_rad = real(xi_width/180.d0 * PI,kind=CUSTOM_REAL)
  eta_width_rad = real(eta_width/180.d0 * PI,kind=CUSTOM_REAL)

  ! gets process/slice number
  xi1 = xi / xi_width_rad * 2.0
  eta1 = eta / eta_width_rad * 2.0

  proc_xi = floor((xi1 + 1.0)/2.0 * nproc_xi)
  proc_eta = floor((eta1 + 1.0)/2.0 * nproc_xi)

  slice_proc = nproc_xi * nproc_xi * (chunk - 1) + nproc_xi * proc_eta + proc_xi

  !print *,"xi1  = ",xi1," xi_width = ",xi_width
  !print *,"eta1 = ",eta1," eta_width = ",eta_width
  !print *,"proc_xi = ",proc_xi
  !print *,"proc_eta = ",proc_eta

  end subroutine get_process_slice_number

