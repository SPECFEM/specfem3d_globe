!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
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
!Vala Hjorleifsdottir, vala@geofisica.unam.mx
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

  use meshfem3D_models_par
  use shared_parameters

  implicit none

!-------------------------------------------------------
! USER parameters

  ! initial position
  double precision,parameter :: COLAT_0 = 1.0
  double precision,parameter :: LON_0   = 1.0

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
  ! for ellipticity
  !integer :: nspl
  !double precision :: rspl(NR),espl(NR),espl2(NR)
  ! proc numbers for MPI
  integer :: sizeprocs

  ! arguments
  integer :: count
  character(len=MAX_STRING_LEN) :: arg

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

  character(len=MAX_STRING_LEN) :: outfile
  character(len=7) :: str_info

! ************** PROGRAM STARTS HERE **************

  ! initializes MPI to be able to call meshfem routines below
  ! starts mpi
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  ! checks number of processes
  ! note: must run as a single process with: mpirun -np 1 ..
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
  initial_colat =  COLAT_0! 90.d0 - start_lon
  initial_lon = LON_0

  ! loop range
  icolat_start = COLAT_istart
  icolat_end = COLAT_iend
  ilon_start = LON_istart
  ilon_end = LON_iend

  ! checks program arguments
  count = command_argument_count()
  if (count > 0) then
    if (count /= 2 .and. count /= 4) call usage()
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
  endif
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
  if (ONE_CRUST) then
    NUMBER_OF_MESH_LAYERS = MAX_NUMBER_OF_MESH_LAYERS - 1
  else
    NUMBER_OF_MESH_LAYERS = MAX_NUMBER_OF_MESH_LAYERS
  endif
!> end of part from create_regions_mesh -> initialize_layers()

  ! loop over all theta (colatitude) and phis (longitude), every two degrees
  do i = icolat_start,icolat_end
    do j = ilon_start,ilon_end
      ! target location
      theta_degrees = initial_colat + i*delta_lat ! colatitude [0,180]
      phi_degrees   = initial_lon   + j*delta_lon ! longitude [0,360]

      ! checks limits
      if (theta_degrees < 0.0) stop 'Error invalid colatitude < 0'
      if (theta_degrees > 180.0) stop 'Error invalid colatitude > 180'
      if (phi_degrees < 0.0) phi_degrees = phi_degrees + 360.d0
      if (phi_degrees > 360.0) phi_degrees = phi_degrees - 360.d0

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
      write(57,'(a)') '# model         : '//trim(MODEL)
      write(57,'(a,F10.4)') '# point location: colatitude [degree] = ',sngl(theta_degrees)
      write(57,'(a,F10.4)') '# point location: longitude  [degree] = ',sngl(phi_degrees)
      if (CRUSTAL) then
        write(57,'(a,F10.4)') '# moho depth    : moho [km] = ',sngl(moho*R_EARTH_KM)
      else
        write(57,'(a,F10.4)') '# moho 1D depth : moho [km] =',sngl((1.d0-RMOHO/R_EARTH)*R_EARTH_KM)
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
      if (OCEANS .or. ((.not. CRUSTAL) .and. (ROCEAN < R_EARTH))) then
        if (OCEANS) then
          write(57,'(a)') '# oceans'
        else
          ! 1D crustal models
          if (TOPOGRAPHY) then
            write(57,'(a)') '# oceans'
          else
            write(57,'(a)') '# 1D reference oceans'
          endif
        endif
      else
        if (TOPOGRAPHY) then
          write(57,'(a)') '# oceans'
        else
          write(57,'(a)') '# no oceans'
        endif
      endif
      write(57,'(a)') '# no ellipticity'
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
          if (rmin == RMOHO_FICTITIOUS_IN_MESHER/R_EARTH) rmin = 1.0d0 - moho
          if (rmax == RMOHO_FICTITIOUS_IN_MESHER/R_EARTH) rmax = 1.0d0 - moho
          !print *,'rmin == moho at line ',iline
        endif

        if (rmin == rmax_last) then  !!!! this means that we have just jumped between layers
          ! depth increment
          ! write values every 10 km in the deep earth and every 1 km in the shallow earth
          if (rmin > (R_EARTH_KM-DELTA_HIRES_DEPTH)/R_EARTH_KM) then
            delta = 1.0d0   ! in km
          else
            delta = 10.0d0  ! in km
          endif

          ! normalization
          delta = delta/R_EARTH_KM

          ! sets maximum radius without ocean for 1D models
          if (((.not. CRUSTAL) .and. (ROCEAN < R_EARTH)) .and. (.not. TOPOGRAPHY)) then
            ! stops at ocean depth and adds last ocean layers explicitly
            if (rmax == 1.0d0) rmax = ROCEAN/R_EARTH
          endif

          ! backup to detect jump between layers
          rmax_last = rmax

          ! number of iterations in increments of delta between rmin and rmax
          nit = floor((rmax - rmin)/delta) + 1

          do idep = 1,nit+1
            ! line counters
            ! inner core boundary
            if (rmin == RICB/R_EARTH .and. idep == 1) iline_icb = iline
            ! core mantle boundary
            if (rmin == RCMB/R_EARTH .and. idep == 1) iline_cmb = iline
            ! moho
            if (CRUSTAL) then
              ! uses 3D crustal model (e.g. Crust2.0)
              if (rmin == (1.0d0 - moho) .and. idep == 1) then
                iline_moho = iline
              endif
            else
              ! 1D crust from reference model
              if (rmin == RMOHO/R_EARTH .and. idep == 1) iline_moho = iline
            endif

            ! radius
            r = rmin + (idep-1)*delta

            ! make sure we are within the right shell in PREM to honor discontinuities
            ! use small geometrical tolerance
            r_prem = r
            if (r <= rmin*1.000001d0) r_prem = rmin*1.000001d0
            if (r >= rmax*0.999999d0) r_prem = rmax*0.999999d0

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

            ! note: ELLIPTICITY will be ignored for now...

            ! scale values read from routines back to true values
            scaleval = dsqrt(PI*GRAV*RHOAV)
            rho = rho*RHOAV/1000.0d0              ! [kg/m3]
            vpv  = vpv*R_EARTH*scaleval/1000.0d0  ! [m/s]
            vph  = vph*R_EARTH*scaleval/1000.0d0  ! [m/s]
            vsv  = vsv*R_EARTH*scaleval/1000.0d0  ! [m/s]
            vsh  = vsh*R_EARTH*scaleval/1000.0d0  ! [m/s]

            ! finally write the values obtained at the given depth to file
            str_info = ''
            if (iline == iline_icb) str_info = ' # ICB'
            if (iline == iline_cmb) str_info = ' # CMB'
            if (iline == iline_moho) str_info = ' # moho'

            write(57,'(F8.0,7F9.2,F9.5,a)') &
              sngl(r_prem*R_EARTH),sngl(rho*1000.d0),sngl(vpv*1000.d0),sngl(vsv*1000.d0), &
              sngl(Qkappa),sngl(Qmu),sngl(vph*1000.d0),sngl(vsh*1000.d0),sngl(eta_aniso),trim(str_info)

            ! line counter
            iline = iline + 1

            ! debug
            !write(*,'(i3,11F10.4)') &
            ! iline,sngl(rmin*R_EARTH_KM),sngl(rmax*R_EARTH_KM),sngl(r_prem*R_EARTH_KM),sngl(r*R_EARTH_KM), &
            ! sngl(vpv),sngl(vph),sngl(vsv),sngl(vsh),sngl(rho),sngl(eta_aniso),sngl(Qmu)

          enddo !idep
        endif !rmin == rmax_last
      enddo !ilayer

      ! adds ocean layer profile
      call write_profile_ocean(r_prem,elevation,iline,iline_ocean)

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
    write(IMAIN,*) 'see profiles written in directory: ',trim(OUTPUT_FILES)
    write(IMAIN,*)
    write(IMAIN,*) 'End of xwrite_profile'
    write(IMAIN,*)
    call flush_IMAIN()
    close(IMAIN)
  endif

  call finalize_mpi()

  end program xwrite_profile

!
!-------------------------------------------------------------------------------------
!

  subroutine usage()

  use meshfem3D_models_par, only: myrank

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

  use meshfem3D_models_par
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
    if (ISOTROPIC_3D_MANTLE) then
      write(IMAIN,*) '  incorporating 3-D lateral variations'
    else
      write(IMAIN,*) '  no 3-D lateral variations'
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
    write(IMAIN,*)
    write(IMAIN,*)
  endif

  if (myrank == 0) print *,'reading model...'

  ! loads 3D models
  call meshfem3D_models_broadcast(NSPEC_REGIONS,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD, &
                                  R80,R220,R670,RCMB,RICB, &
                                  LOCAL_PATH)

  ! user output
  if (myrank == 0) print *

!> end part from meshfem3D

  end subroutine write_profile_setup

!
!-------------------------------------------------------------------------------------
!


  subroutine write_profile_moho_depth(theta,phi,moho)

  use meshfem3D_models_par
  use shared_parameters

  implicit none

  double precision,intent(in) :: theta,phi
  double precision,intent(out) :: moho

  ! local parameters
  double precision :: r,xmesh,ymesh,zmesh
  double precision :: dvp,rho,vpv,vph,vsv,vsh,eta_aniso
  double precision :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  ! moho depth
  ! read crustal models and topo models as they are needed to modify the depths of the discontinuities
  if (CRUSTAL) then
    ! convert from rthetaphi to xyz, with r=-7km
    r = 1.0d0 - 7.0d0/R_EARTH_KM
    call rthetaphi_2_xyz_dble(xmesh,ymesh,zmesh,r,theta,phi)

    !if (myrank == 0) print *, '  xmesh,ymesh,zmesh,r,theta,phi = ',xmesh,ymesh,zmesh,r,90.0-theta*180./PI,phi*180./PI

    ! gets moho depth from crustal model
    call meshfem3D_models_get3Dcrust_val(IREGION_CRUST_MANTLE,xmesh,ymesh,zmesh,r, &
                                         vpv,vph,vsv,vsh,rho,eta_aniso,dvp, &
                                         c11,c12,c13,c14,c15,c16,c22,c23,c24,c25, &
                                         c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66, &
                                         .true.,moho)

    !if (myrank == 0) print *, '  crustal values:',vpv,vph,vsv,vsh,rho,eta_aniso,dvp
    if (myrank == 0) print *, '  moho depth [km]:',sngl(moho*R_EARTH_KM) !, &
                              !'moho radius (in km):',sngl((1.0d0-moho)*R_EARTH_KM)
  else
    moho = 0.d0
    ! 1D crust from reference model
    if (myrank == 0) print *, '  moho 1D depth [km]:',sngl((1.d0-RMOHO/R_EARTH)*R_EARTH_KM)
  endif

  end subroutine write_profile_moho_depth

!
!-------------------------------------------------------------------------------------
!

  subroutine write_profile_elevation(theta,phi,elevation)

  use meshfem3D_models_par
  use shared_parameters

  implicit none

  double precision,intent(in) :: theta,phi
  double precision,intent(out) :: elevation

  ! local parameters
  double precision :: lat,lon

  ! topography elevation
  if (TOPOGRAPHY .or. OCEANS) then
    if (TOPOGRAPHY) then
      lat = (PI/2.0d0-theta)*180.0d0/PI     ! range [-90,90]
      lon = phi*180.0d0/PI
      if (lon > 180.0d0) lon = lon-360.0d0  ! range [-180,180]

      !if (myrank == 0) print *,'  get_topo_bathy(lat,lon,elevation,ibathy_topo',lat,lon,elevation

      call get_topo_bathy(lat,lon,elevation,ibathy_topo)

      !if (myrank == 0) print *,'  get_topo_bathy(lat,lon,elevation,ibathy_topo',lat,lon,elevation
      if (myrank == 0) print *,'  elevation [m]:',sngl(elevation) !, &
                               !'surface radius (in km):',sngl((1.0d0 + elevation/R_EARTH)*R_EARTH_KM)

    else
      ! oceans selected only, without topography
      if ((.not. CRUSTAL) .and. (ROCEAN < R_EARTH)) then
        ! if 1D Earth, use oceans of constant thickness everywhere (see create_mass_matrices_ocean_load())
        elevation = - THICKNESS_OCEANS_PREM * R_EARTH ! in m
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

  use meshfem3D_models_par
  use shared_parameters

  implicit none

  double precision,intent(in) :: r,r_prem,theta,phi
  integer,intent(in) :: iregion_code,idoubling
  double precision,intent(in) :: rmin,rmax
  double precision,intent(out) :: rho,vpv,vph,vsv,vsh,eta_aniso,Qmu,Qkappa

  ! local parameters
  double precision :: xmesh,ymesh,zmesh
  double precision :: dvp,moho
  double precision :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66
  ! Attenuation values
  double precision, dimension(N_SLS) :: tau_s, tau_e
  double precision :: T_c_source
  ! Parameter used to decide whether this element is in the crust or not
  logical:: elem_in_crust,elem_in_mantle

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
  dvp = 0.d0

  ! do not force the element/point to be in the crust
  elem_in_crust = .false.

  ! convert from rthetaphi to xyz to use in function calls.
  call rthetaphi_2_xyz_dble(xmesh,ymesh,zmesh,r_prem,theta,phi)

! < start GET_MODEL
  ! checks r_prem,rmin/rmax and assigned idoubling
  call get_model_check_idoubling(r_prem,xmesh,ymesh,zmesh,rmin,rmax,idoubling, &
                      RICB,RCMB,RTOPDDOUBLEPRIME, &
                      R220,R670)

  ! gets reference model values: rho,vpv,vph,vsv,vsh and eta_aniso
  call meshfem3D_models_get1D_val(iregion_code,idoubling, &
                        r_prem,rho,vpv,vph,vsv,vsh,eta_aniso, &
                        Qkappa,Qmu,RICB,RCMB, &
                        RTOPDDOUBLEPRIME,R80,R120,R220,R400,R600,R670,R771, &
                        RMOHO,RMIDDLE_CRUST,ROCEAN)

  ! gets the 3-D model parameters for the mantle
  call meshfem3D_models_get3Dmntl_val(iregion_code,r_prem,rho,dvp, &
                        vpv,vph,vsv,vsh,eta_aniso, &
                        RCMB,R670,RMOHO, &
                        xmesh,ymesh,zmesh,r, &
                        c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                        c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

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
      call meshfem3D_models_get3Dcrust_val(iregion_code,xmesh,ymesh,zmesh,r_prem, &
                        vpv,vph,vsv,vsh,rho,eta_aniso,dvp, &
                        c11,c12,c13,c14,c15,c16,c22,c23,c24,c25, &
                        c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66, &
                        elem_in_crust,moho)

      !print *,'crustal: ',r_prem,R_DEEPEST_CRUST,elem_in_crust,vpv,vph,vsv,vsh,eta_aniso
    endif
  endif

  !!! VH  commented out following two lines from get_model
  !! overwrites with tomographic model values (from iteration step) here, given at all GLL points
  !call meshfem3D_models_impose_val(vpv,vph,vsv,vsh,rho,dvp,eta_aniso,iregion_code,ispec,i,j,k)

  ! checks vpv: if close to zero then there is probably an error
  if (vpv < TINYVAL) then
    print *,'error vpv: ',vpv,vph,vsv,vsh,rho
    print *,'radius:',r*R_EARTH_KM
    call exit_mpi(myrank,'Error get model values in write_profile')
  endif

  !> Hejun
  ! New Attenuation assignment
  ! Define 3D and 1D Attenuation after moho stretch
  ! and before TOPOGRAPHY/ELLIPCITY
  !
  !note:  only Qmu attenuation considered, Qkappa attenuation not used so far...
  if (ATTENUATION) then
    call meshfem3D_models_getatten_val(idoubling,xmesh,ymesh,zmesh,r_prem, &
                                       tau_e,tau_s,T_c_source, &
                                       moho,Qmu,Qkappa,elem_in_crust)
  endif
!> end GET_MODEL

  end subroutine write_profile_model_values

!
!-------------------------------------------------------------------------------------
!


  subroutine write_profile_add_topography(r_prem,idoubling,elevation)

  use meshfem3D_models_par
  use shared_parameters

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
      gamma = (r_prem - R220/R_EARTH) / (R_UNIT_SPHERE - R220/R_EARTH)
      if (gamma < -0.02 .or. gamma > 1.02) print *, 'incorrect value of gamma for topograpy'

      ! stretching due to elevation
      !print *,'rprem before: ',r_prem*R_EARTH
      r_prem = r_prem*(ONE + gamma * (elevation/R_EARTH) /r_prem)
      !print *,'r_prem after: ',r_prem*R_EARTH
    endif
  endif
!> end add_topography

  end subroutine write_profile_add_topography

!
!-------------------------------------------------------------------------------------
!


  subroutine write_profile_ocean(r_prem,elevation,iline,iline_ocean)

  use meshfem3D_models_par
  use shared_parameters

  implicit none

  double precision,intent(in) :: r_prem,elevation
  integer,intent(inout) :: iline,iline_ocean

  ! local parameters
  ! depth increment
  double precision :: delta

  ! ocean
  integer :: ilayers_ocean,nlayers_ocean
  double precision :: r_ocean

! < ocean
  ! This part adds the ocean to profile where needed
  if ( (TOPOGRAPHY .or. OCEANS) .or. ((.not. CRUSTAL) .and. (ROCEAN < R_EARTH)) ) then

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
        if ((.not. CRUSTAL) .and. (ROCEAN < R_EARTH)) then
          ! 1D model without topography, uses default ocean layer (R_EARTH - ROCEAN)
          nlayers_ocean = floor((R_EARTH - ROCEAN)/delta)
        else
          ! 3D crustal model, without topography
          ! no oceans
          nlayers_ocean = -1
        endif
      endif
    else
      ! case for 1D (reference crust) models which include an ocean: ((.not. CRUSTAL) .and. (ROCEAN < R_EARTH)) then
      if (TOPOGRAPHY) then
        if (elevation < - MINIMUM_THICKNESS_3D_OCEANS) then
          ! uses actual bathymetry value and fills up ocean to surface
          nlayers_ocean = floor(-elevation/delta) ! this means ocean depth must be bigger than 500m (delta) for output
        else
          ! elevation above minimum ocean depth
          nlayers_ocean = -1
        endif
      else
        if ((.not. CRUSTAL) .and. (ROCEAN < R_EARTH)) then
          ! 1D models without topography, takes reference ocean depth
          nlayers_ocean = floor((R_EARTH - ROCEAN)/delta)
        else
          ! 3D crustal models without topography will have no ocean
          nlayers_ocean = -1
        endif
      endif
    endif

    ! adds ocean layer
    if (nlayers_ocean > 0) then
      ! ocean line
      iline_ocean = iline

      ! points within ocean
      do ilayers_ocean = 0,nlayers_ocean
        ! radius
        r_ocean = r_prem + ilayers_ocean* delta/R_EARTH

        ! ocean properties (salt water parameters from PREM)
        if (ilayers_ocean == 0) then
          ! line with section info
          write(57,'(F8.0,7F9.2,F9.5,a)') sngl(r_ocean*R_EARTH),1020.0,1450.,0.0,57822.5,0.0,1450.0,0.0,1.0,' # ocean'
        else
          write(57,'(F8.0,7F9.2,F9.5)') sngl(r_ocean*R_EARTH),1020.0,1450.,0.0,57822.5,0.0,1450.0,0.0,1.0
        endif
        ! line counter
        iline = iline + 1
      enddo
      ! at surface
      if (r_ocean < 1.d0) then
        ! last line exactly at earth surface
        write(57,'(F8.0,7F9.2,F9.5)') sngl(1.0d0*R_EARTH),1020.0,1450.,0.0,57822.5,0.0,1450.,0.0,1.0
        ! line counter
        iline = iline + 1
      endif
    endif
  endif
!> end adding ocean

  end subroutine write_profile_ocean
