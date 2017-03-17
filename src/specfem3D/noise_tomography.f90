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

! subroutines for NOISE TOMOGRAPHY


  subroutine noise_distribution_direction(r_in,theta_in,phi_in, &
                                          normal_x_noise_out,normal_y_noise_out,normal_z_noise_out, &
                                          mask_noise_out)

! characterizes noise statistics
!
! for a given point (xcoord,ycoord,zcoord), specify the noise direction "normal_x/y/z_noise"
!     and noise distribution "mask_noise"
! USERS need to modify this subroutine for their own noise characteristics

  use constants

  implicit none

  ! input parameters
  real(kind=CUSTOM_REAL) :: r_in,theta_in,phi_in
  ! output parameters
  real(kind=CUSTOM_REAL) :: normal_x_noise_out,normal_y_noise_out,normal_z_noise_out,mask_noise_out
  ! local parameters
  real(kind=CUSTOM_REAL) :: xcoord,ycoord,zcoord


  ! convert r theta phi back to x y z
  call rthetaphi_2_xyz(xcoord,ycoord,zcoord,r_in,theta_in,phi_in)
  ! NOTE that all coordinates are non-dimensionalized in GLOBAL package!
  ! USERS are free to choose which set to use,
  ! either "r theta phi" (xcoord_in,ycoord_in,zcoord_in)
  ! or     "x y z"       (xcoord,ycoord,zcoord)

  !*****************************************************************************************************************
  !******************************** change your noise characteristics below ****************************************
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! noise direction
  ! here, the noise is assumed to be vertical
  normal_x_noise_out = xcoord / sqrt(xcoord**2 + ycoord**2 + zcoord**2)
  normal_y_noise_out = ycoord / sqrt(xcoord**2 + ycoord**2 + zcoord**2)
  normal_z_noise_out = zcoord / sqrt(xcoord**2 + ycoord**2 + zcoord**2)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  noise distribution
  ! here, the noise is assumed to be uniform
  mask_noise_out = 1.0
  !******************************** change your noise characteristics above ****************************************
  !*****************************************************************************************************************

  end subroutine noise_distribution_direction

!
!-------------------------------------------------------------------------------------------------
!

  subroutine check_parameters_noise()

! checks for consistency of the parameters

  use specfem_par
  use specfem_par_crustmantle, only: NSPEC_TOP
  use specfem_par_noise

  implicit none

  ! local parameters
  integer :: ier
  character(len=MAX_STRING_LEN) :: filename

  ! checks if anything to do
  if (NOISE_TOMOGRAPHY == 0) return

  ! info file output
  if (myrank == 0) then
    filename = trim(OUTPUT_FILES)//'/NOISE_SIMULATION.txt'
    open(unit=IOUT_NOISE,file=trim(filename),status='unknown',action='write',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening output file NOISE_SIMULATION.txt')

    write(IOUT_NOISE,*) '*******************************************************************************'
    write(IOUT_NOISE,*) '*******************************************************************************'
    write(IOUT_NOISE,*) 'WARNING!!!!!!!!!!!!'
    write(IOUT_NOISE,*) 'You are running simulations using NOISE TOMOGRAPHY techniques.'
    write(IOUT_NOISE,*) 'Please make sure you understand the procedures before you have a try.'
    write(IOUT_NOISE,*)
    write(IOUT_NOISE,*) 'Displacements everywhere at the free surface are saved every timestep,'
    write(IOUT_NOISE,*) 'so make sure that LOCAL_TMP_PATH in DATA/Par_file is not global.'
    write(IOUT_NOISE,*) 'Otherwise the disk storage may be a serious issue, as is the speed of I/O.'
    write(IOUT_NOISE,*)
    write(IOUT_NOISE,*) 'Also note that NO earthquakes are included,'
    write(IOUT_NOISE,*) 'i.e., the moment tensor in CMTSOLUTION will be ignored.'
    write(IOUT_NOISE,*)
    write(IOUT_NOISE,*) 'If you just want a regular EARTHQUAKE simulation,'
    write(IOUT_NOISE,*) 'set NOISE_TOMOGRAPHY = 0 in DATA/Par_file'
    write(IOUT_NOISE,*) '*******************************************************************************'
    write(IOUT_NOISE,*) '*******************************************************************************'
    close(IOUT_NOISE)
  endif

  ! checks parameters
  if (NUMBER_OF_RUNS /= 1 .or. NUMBER_OF_THIS_RUN /= 1) &
    call exit_mpi(myrank,'NUMBER_OF_RUNS and NUMBER_OF_THIS_RUN must be 1 for NOISE TOMOGRAPHY! check DATA/Par_file')
  if (ROTATE_SEISMOGRAMS_RT) &
    call exit_mpi(myrank,'Do NOT rotate seismograms in the code, change ROTATE_SEISMOGRAMS_RT in DATA/Par_file')
  if (SAVE_ALL_SEISMOS_IN_ONE_FILE .or. USE_BINARY_FOR_LARGE_FILE) &
    call exit_mpi(myrank,'Please set SAVE_ALL_SEISMOS_IN_ONE_FILE and USE_BINARY_FOR_LARGE_FILE to be .false.')

  ! well, experimental, but free to try out...
  !if (UNDO_ATTENUATION .and. NOISE_TOMOGRAPHY > 0 ) &
  !  call exit_mpi(myrank,'UNDO_ATTENUATION support not implemented yet for noise simulations')

  ! checks noise simulation setup
  select case (NOISE_TOMOGRAPHY)
  case (1)
    ! forward noise source simulation
    if (SIMULATION_TYPE /= 1) &
      call exit_mpi(myrank,'NOISE_TOMOGRAPHY=1 requires SIMULATION_TYPE=1! check DATA/Par_file')
  case (2)
    ! forward ensemble wavefield simulation
    if (SIMULATION_TYPE /= 1) &
      call exit_mpi(myrank,'NOISE_TOMOGRAPHY=2 requires SIMULATION_TYPE=1! check DATA/Par_file')
    if (.not. SAVE_FORWARD) &
      call exit_mpi(myrank,'NOISE_TOMOGRAPHY=2 requires SAVE_FORWARD=.true.! check DATA/Par_file')
  case (3)
    ! adjoint ensemble kernel simulation
    if (SIMULATION_TYPE /= 3) &
      call exit_mpi(myrank,'NOISE_TOMOGRAPHY=3 requires SIMULATION_TYPE=3! check DATA/Par_file')
    if (SAVE_FORWARD) &
      call exit_mpi(myrank,'NOISE_TOMOGRAPHY=3 requires SAVE_FORWARD=.false.! check DATA/Par_file')
  case default
    call exit_MPI(myrank,'Error invalid NOISE_TOMOGRAPHY value for noise simulation setup! check DATA/Par_file')
  end select

  ! checks that number of spectral elements at surface is set (from read_mesh_databases_CM() routine)
  if (NSPEC_TOP /= NSPEC2D_TOP(IREGION_CRUST_MANTLE)) &
    call exit_MPI(myrank,'Error invalid number of NSPEC_TOP for noise simulation')

  ! save/read the surface movie using the same c routine as we do for absorbing boundaries (file ID is 9)
  ! size of single record
  reclen_noise = CUSTOM_REAL * NDIM * NGLLX * NGLLY * NSPEC_TOP

  ! check integer size limit: size of reclen_noise must fit onto an 4-byte integer
  if (NSPEC_TOP > 2147483646 / (CUSTOM_REAL * NGLLX * NGLLY * NDIM)) then
    print *,'size of noise surface movie array needed exceeds integer 4-byte limit: ', &
      dble(CUSTOM_REAL * NGLLX * NGLLY * NDIM) * dble(NSPEC_TOP),reclen_noise
    print *,'  ',CUSTOM_REAL, NDIM, NGLLX * NGLLY, NSPEC_TOP
    print *,'bit size fortran: ',bit_size(NSPEC_TOP)
    call exit_MPI(myrank,"Error NSPEC_TOP leads to noise surface array exceeding integer limit (2 GB)")
  endif

  end subroutine check_parameters_noise

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_parameters_noise()

! reads noise parameters

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_noise

  implicit none

  ! local parameters
  integer :: ipoin, ispec2D, ispec, i, j, k, iglob, ier
  integer(kind=8) :: filesize
  real(kind=CUSTOM_REAL) :: normal_x_noise_out,normal_y_noise_out,normal_z_noise_out,mask_noise_out
  real(kind=CUSTOM_REAL) :: r,theta,phi
  character(len=MAX_STRING_LEN) :: filename

  ! read master receiver ID -- the ID in DATA/STATIONS
  filename = trim(OUTPUT_FILES)//'/..//NOISE_TOMOGRAPHY/irec_master_noise'
  open(unit=IIN_NOISE,file=trim(filename),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    call exit_MPI(myrank, 'file '//trim(filename)//' does NOT exist! This file contains the ID of the master receiver')
  endif

  read(IIN_NOISE,*,iostat=ier) irec_master_noise
  if (ier /= 0) then
    call exit_MPI(myrank, 'Unable to read the ID of the master receiver from '//trim(filename))
  endif
  close(IIN_NOISE)

  ! user output
  if (myrank == 0) then
    filename = trim(OUTPUT_FILES)//'/irec_master_noise'
    open(unit=IOUT_NOISE,file=trim(filename),status='unknown',action='write',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening output file irec_master_noise')
    write(IOUT_NOISE,*) 'The master receiver is: (RECEIVER ID)', irec_master_noise
    close(IOUT_NOISE)
  endif

  ! sets number of noise master sources located in this slice (only 1 slice possible to hold master station)
  ! (only used and needed for 1. step of noise simulations)
  nsources_local_noise = 0
  if (NOISE_TOMOGRAPHY == 1) then
    ! checks master irec
    if (irec_master_noise < 1 .or. irec_master_noise > nrec) then
      call exit_MPI(myrank,'Error noise tomography: irec_master_noise is not in range of given number of receivers')
    endif

    ! increases to one for slice holding master stations
    if (myrank == islice_selected_rec(irec_master_noise)) nsources_local_noise = nsources_local_noise + 1

    ! compute source arrays for "ensemble forward source", which is source of "ensemble forward wavefield"
    if (nsources_local_noise > 0 .or. myrank == 0) then ! myrank == 0 is used for output only
      call compute_arrays_source_noise(xi_receiver(irec_master_noise), &
                                       eta_receiver(irec_master_noise), &
                                       gamma_receiver(irec_master_noise), &
                                       nu(:,:,irec_master_noise),noise_sourcearray, xigll,yigll,zigll,NSTEP)
    endif
  endif


  ! sets up noise distribution and noise direction
  !
  ! loops over all surface points
  ! puts noise distrubution and direction onto the surface points
  ipoin = 0
  do ispec2D = 1, NSPEC_TOP ! NSPEC2D_TOP(IREGION_CRUST_MANTLE)
    ispec = ibelm_top_crust_mantle(ispec2D)

    ! loop on all the points at the surface of the element
    k = NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        ipoin = ipoin + 1
        iglob = ibool_crust_mantle(i,j,k,ispec)

        ! point location
        r = rstore_crust_mantle(1,iglob)
        theta = rstore_crust_mantle(2,iglob)
        phi = rstore_crust_mantle(3,iglob)

        ! this subroutine must be modified by USERS
        call noise_distribution_direction(r,theta,phi, &
                                          normal_x_noise_out,normal_y_noise_out,normal_z_noise_out, &
                                          mask_noise_out)

        ! stores normal components
        normal_x_noise(ipoin) = normal_x_noise_out
        normal_y_noise(ipoin) = normal_y_noise_out
        normal_z_noise(ipoin) = normal_z_noise_out

        ! stores point mask
        mask_noise(ipoin) = mask_noise_out
      enddo
    enddo

  enddo
  ! checks
  if (ipoin /= num_noise_surface_points) call exit_MPI(myrank,'Error invalid number of surface points for noise distribution')

  ! saves mask_noise for check, a file called "mask_noise.bin" is saved in "./OUTPUT_FILES/"
  call save_mask_noise()

  ! opens noise surface movie files for file I/O
  !
  ! total file size
  filesize = reclen_noise
  filesize = filesize * NSTEP

  ! noise surface array stored by each process
  write(filename,"('/proc',i6.6,'_noise_surface_movie.bin')") myrank

  if (NOISE_TOMOGRAPHY == 1) call open_file_abs_w(9,trim(LOCAL_TMP_PATH)//trim(filename), &
                                                  len_trim(trim(LOCAL_TMP_PATH)//trim(filename)), filesize)
  if (NOISE_TOMOGRAPHY == 2) call open_file_abs_r(9,trim(LOCAL_TMP_PATH)//trim(filename), &
                                                  len_trim(trim(LOCAL_TMP_PATH)//trim(filename)), filesize)
  if (NOISE_TOMOGRAPHY == 3) call open_file_abs_r(9,trim(LOCAL_TMP_PATH)//trim(filename), &
                                                  len_trim(trim(LOCAL_TMP_PATH)//trim(filename)), filesize)

  end subroutine read_parameters_noise

!
!-------------------------------------------------------------------------------------------------
!

  subroutine save_mask_noise()

! saves mask_noise array for check, a file called "mask_noise.bin" is saved in "./OUTPUT_FILES/"

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_noise

  implicit none

  ! local parameters
  integer :: ipoin, ispec2D, ispec, i, j, k, iglob, ier
  character(len=MAX_STRING_LEN) :: filename
  real(kind=CUSTOM_REAL), dimension(num_noise_surface_points) :: &
      val_x,val_y,val_z,val_ux,val_uy,val_uz
  real(kind=CUSTOM_REAL), dimension(num_noise_surface_points,0:NPROCTOT_VAL-1) :: &
      val_x_all,val_y_all,val_z_all,val_ux_all,val_uy_all,val_uz_all

  ! sets up temporary arrays for this slice
  ipoin = 0
  do ispec2D = 1, NSPEC_TOP ! NSPEC2D_TOP(IREGION_CRUST_MANTLE)
    ispec = ibelm_top_crust_mantle(ispec2D)
    k = NGLLZ
    ! loop on all the points inside the element
    do j = 1,NGLLY
      do i = 1,NGLLX
        ipoin = ipoin + 1
        iglob = ibool_crust_mantle(i,j,k,ispec)
        ! stores position
        val_x(ipoin) = rstore_crust_mantle(1,iglob) ! r
        val_y(ipoin) = rstore_crust_mantle(2,iglob) ! theta
        val_z(ipoin) = rstore_crust_mantle(3,iglob) ! phi
        ! stores mask
        val_ux(ipoin) = mask_noise(ipoin)
        val_uy(ipoin) = mask_noise(ipoin)
        val_uz(ipoin) = mask_noise(ipoin)
      enddo
    enddo
  enddo
  ! checks
  if (ipoin /= num_noise_surface_points) call exit_MPI(myrank,'Error invalid number of surface points for noise mask')

  ! gather info on master proc
  call gather_all_cr(val_x,num_noise_surface_points,val_x_all,num_noise_surface_points,NPROCTOT_VAL)
  call gather_all_cr(val_y,num_noise_surface_points,val_y_all,num_noise_surface_points,NPROCTOT_VAL)
  call gather_all_cr(val_z,num_noise_surface_points,val_z_all,num_noise_surface_points,NPROCTOT_VAL)

  call gather_all_cr(val_ux,num_noise_surface_points,val_ux_all,num_noise_surface_points,NPROCTOT_VAL)
  call gather_all_cr(val_uy,num_noise_surface_points,val_uy_all,num_noise_surface_points,NPROCTOT_VAL)
  call gather_all_cr(val_uz,num_noise_surface_points,val_uz_all,num_noise_surface_points,NPROCTOT_VAL)

  ! saves mask_noise data to disk in home directory
  !
  ! this file can be viewed the same way as surface movie data (xcreate_movie_AVS_DX)
  ! create_movie_AVS_DX.f90 needs to be modified in order to do that,
  ! i.e., instead of showing the normal component, change it to either x, y or z component, or the norm.
  if (myrank == 0) then
    ! a file called "mask_noise.bin" is saved in "./OUTPUT_FILES/"
    filename = trim(OUTPUT_FILES)//'/mask_noise.bin'

    open(unit=IOUT_NOISE,file=trim(filename),status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening output file mask_noise')

    write(IOUT_NOISE) val_x_all ! r
    write(IOUT_NOISE) val_y_all ! theta
    write(IOUT_NOISE) val_z_all ! phi
    write(IOUT_NOISE) val_ux_all
    write(IOUT_NOISE) val_uy_all
    write(IOUT_NOISE) val_uz_all

    close(IOUT_NOISE)
  endif

  end subroutine save_mask_noise

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_arrays_source_noise(xi_noise,eta_noise,gamma_noise,nu_single,noise_sourcearray, &
                                         xigll,yigll,zigll,NSTEP)

! reads and constructs the "source" (source time function based upon noise spectrum)
! for "ensemble forward source"

  use constants_solver
  use shared_parameters, only: OUTPUT_FILES,DT

  implicit none

  ! input parameters
  integer :: NSTEP
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll
  double precision, dimension(NDIM,NDIM) :: nu_single  ! rotation matrix at the master receiver
  ! output parameters
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ,NSTEP) :: noise_sourcearray
  ! local parameters
  integer itime, i, j, k, ier
  real(kind=CUSTOM_REAL) :: junk
  real(kind=CUSTOM_REAL), dimension(NSTEP) :: noise_src
  real(kind=CUSTOM_REAL), dimension(NDIM,NSTEP) :: noise_src_u
  double precision, dimension(NDIM) :: nu_master       ! component direction chosen at the master receiver
  double precision :: xi_noise, eta_noise, gamma_noise ! master receiver location
  double precision, dimension(NGLLX) :: hxir, hpxir
  double precision, dimension(NGLLY) :: hetar, hpetar
  double precision, dimension(NGLLZ) :: hgammar, hpgammar
  character(len=MAX_STRING_LEN) :: filename

  ! noise file (source time function)
  filename = trim(OUTPUT_FILES)//'/..//NOISE_TOMOGRAPHY/S_squared'
  open(unit=IIN_NOISE,file=trim(filename),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    call exit_MPI(myrank, 'file '//trim(filename)//' does NOT exist! This file should have been generated using Matlab scripts')
  endif

  noise_src(:) = 0._CUSTOM_REAL
  do itime  = 1,NSTEP
    read(IIN_NOISE,*,iostat=ier) junk, noise_src(itime)
    if (ier /= 0) then
      print *,'Error noise source S_squared file length: NSTEP length required is ',NSTEP,' with time step size ',DT
      call exit_MPI(myrank,'file '//trim(filename)//' has wrong length, please check with your simulation duration')
    endif
  enddo
  close(IIN_NOISE)

  ! master receiver component direction, \nu_master
  filename = trim(OUTPUT_FILES)//'/..//NOISE_TOMOGRAPHY/nu_master'
  open(unit=IIN_NOISE,file=trim(filename),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    call exit_MPI(myrank, &
                  'file '//trim(filename)//' does NOT exist! nu_master is the component direction (NEZ) for master receiver')
  endif

  do itime = 1,3
    read(IIN_NOISE,*,iostat=ier) nu_master(itime)
    if (ier /= 0) then
      print *,'Error noise nu_master file length: number of required components is 3'
      call exit_MPI(myrank, &
                    'file '//trim(filename)//' has wrong length, the vector should have three components (NEZ)')
    endif
  enddo
  close(IIN_NOISE)

  if (myrank == 0) then
     open(unit=IOUT_NOISE,file=trim(OUTPUT_FILES)//'/nu_master',status='unknown',action='write')
     write(IOUT_NOISE,*) 'The direction (NEZ) of selected component of master receiver is', nu_master
     close(IOUT_NOISE)
  endif

  ! rotates to Cartesian
  do itime = 1, NSTEP
    noise_src_u(:,itime) = nu_single(1,:) * noise_src(itime) * nu_master(1) &
                         + nu_single(2,:) * noise_src(itime) * nu_master(2) &
                         + nu_single(3,:) * noise_src(itime) * nu_master(3)
  enddo

  ! receiver interpolators
  call lagrange_any(xi_noise,NGLLX,xigll,hxir,hpxir)
  call lagrange_any(eta_noise,NGLLY,yigll,hetar,hpetar)
  call lagrange_any(gamma_noise,NGLLZ,zigll,hgammar,hpgammar)

  ! adds interpolated source contribution to all GLL points within this element
  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX
        do itime = 1, NSTEP
          noise_sourcearray(:,i,j,k,itime) = hxir(i) * hetar(j) * hgammar(k) * noise_src_u(:,itime)
        enddo
      enddo
    enddo
  enddo

  end subroutine compute_arrays_source_noise

!
!-------------------------------------------------------------------------------------------------
!

  subroutine print_master_distances_noise()

! lists distances to master stations

  use specfem_par, only: nrec,stlat,stlon,station_name,network_name, &
    HUGEVAL,DEGREES_TO_RADIANS,RADIANS_TO_DEGREES,myrank,IMAIN,NOISE_TOMOGRAPHY
  use specfem_par_noise

  implicit none

  ! local parameters
  integer :: irec,i
  double precision :: lat,lon
  double precision :: theta,phi
  double precision :: theta_master,phi_master
  double precision, dimension(nrec) :: epidist
  ! sorting order
  integer, dimension(nrec) :: irec_dist_ordered

  ! checks if anything to do
  if (NOISE_TOMOGRAPHY /= 1) return

  ! lists distances to master station
  ! master station location
  lat = stlat(irec_master_noise)
  lon = stlon(irec_master_noise)

  ! limits longitude to [0.0,360.0]
  if (lon < 0.d0 ) lon = lon + 360.d0
  if (lon > 360.d0 ) lon = lon - 360.d0

  ! converts geographic latitude stlat (degrees) to geocentric colatitude theta (radians)
  call lat_2_geocentric_colat_dble(lat,theta_master)

  phi_master = lon * DEGREES_TO_RADIANS
  call reduce(theta_master,phi_master)

  ! loop on all the stations to locate them in the mesh
  do irec = 1,nrec
    ! station lat/lon in degrees
    lat = stlat(irec)
    lon = stlon(irec)

    ! limits longitude to [0.0,360.0]
    if (lon < 0.d0 ) lon = lon + 360.d0
    if (lon > 360.d0 ) lon = lon - 360.d0

    ! converts geographic latitude stlat (degrees) to geocentric colatitude theta (radians)
    call lat_2_geocentric_colat_dble(lat,theta)

    phi = lon*DEGREES_TO_RADIANS
    call reduce(theta,phi)

    ! computes epicentral distance
    epidist(irec) = acos(cos(theta)*cos(theta_master) + &
                           sin(theta)*sin(theta_master)*cos(phi-phi_master))*RADIANS_TO_DEGREES
  enddo

  ! print some information about stations
  if (myrank == 0) then
    ! sorts stations according to epicentral distances
    ! sorts array
    call heap_sort_distances(nrec,epidist,irec_dist_ordered)

    ! outputs info
    write(IMAIN,*) 'Stations sorted by epicentral distance to noise master station: ', &
                   trim(network_name(irec_master_noise))//'.'//trim(station_name(irec_master_noise))
    do i = 1,nrec
      irec = irec_dist_ordered(i)
      write(IMAIN,'(a,i6,a,a24,a,f12.6,a)') ' Station #',irec,': ', &
        trim(network_name(irec))//'.'//trim(station_name(irec)), &
        '    epicentral distance:  ',sngl(epidist(irec)),' degrees'
    enddo
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine print_master_distances_noise

!
!-------------------------------------------------------------------------------------------------
!

  subroutine noise_add_source_master_rec()

! step 1: calculate the "ensemble forward source"
! add noise spectrum to the location of master receiver
!
! the first step of noise tomography is to use |S(\omega)|^2 as a point force source at one of the receivers.
! hence, instead of a moment tensor 'sourcearrays', a 'noise_sourcearray' for a point force is needed.
! furthermore, the CMTSOLUTION needs to be zero, i.e., no earthquakes.
! now this must be manually set in DATA/CMTSOLUTION, by USERS.
!
! only called for NOISE_TOMOGRAPHY == 1

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_noise

  implicit none

  ! local parameters
  integer :: i,j,k,iglob

  ! checks if anything to do (only 1 slice will contain master station)
  if (nsources_local_noise < 1) return

  ! adds noise source (only if this proc carries the noise)
  if (.not. GPU_MODE) then
    ! on CPU
    ! adds noise source contributions
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool_crust_mantle(i,j,k,ispec_selected_rec(irec_master_noise))
          accel_crust_mantle(:,iglob) = accel_crust_mantle(:,iglob) &
                                        + noise_sourcearray(:,i,j,k,it)
        enddo
      enddo
    enddo
  else
    ! on GPU
    call noise_add_source_master_rec_gpu(Mesh_pointer,it,irec_master_noise,islice_selected_rec)
  endif

  end subroutine noise_add_source_master_rec

!
!-------------------------------------------------------------------------------------------------
!

  subroutine noise_save_surface_movie()

! step 1: calculate the "ensemble forward source"
! save surface movie (displacement) at every time steps, for step 2 & 3.
!
! there are two subroutines --- noise_save_surface_movie_original & noise_save_surface_movie
!    noise_save_surface_movie_original is implemented at first, which creates one file at each time step
!    noise_save_surface_movie is implemented later, which utilizes 'src/shared/binary_c_io.c' for faster I/O,
!                                                   which creates one file for the all time steps
!
! by this modification, the efficiency is greatly improved
! and now, it should be OK to run NOISE_TOMOGRAPHY on a cluster with global storage
!
! only called for NOISE_TOMOGRAPHY == 1

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_noise

  implicit none

  ! local parameters
  integer :: ispec2D,ispec,i,j,k,iglob
  ! i/o buffer
  integer :: it_write_index,it_buffer

  ! get coordinates of surface mesh and surface displacement
  if (.not. GPU_MODE) then
    ! on CPU
    do ispec2D = 1, NSPEC_TOP ! NSPEC2D_TOP(IREGION_CRUST_MANTLE)
      ispec = ibelm_top_crust_mantle(ispec2D)
      k = NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool_crust_mantle(i,j,k,ispec)
          noise_surface_movie(:,i,j,ispec2D) = displ_crust_mantle(:,iglob)
        enddo
      enddo
    enddo
  else
    ! on GPU
    call noise_transfer_surface_to_host(Mesh_pointer,noise_surface_movie)
  endif

  ! save surface motion to disk
  !
  ! original
  !call write_abs(9,noise_surface_movie,reclen_noise,it)
  !
  ! note: due to heavy file i/o, the output of the noise surface movie at every timestep will heavily hit the file servers;
  !       we will be using chunks of multiple time steps to alleviate this.
  !
  !       however, performance of writing will be almost the same, since it only appends data (for binary i/o, not map i/o)
  !       without the need of repositioning the file head. Repositioing will be use for reading in reverse order in
  !       routine noise_read_add_surface_movie() and slow it down considerably.
  !
  ! sets buffer steps and index for file access
  ! (icounter_noise_buffer is a counter between 0 and NT_DUMP_NOISE_BUFFER-1)
  if (icounter_noise_buffer == 0) then
    ! sets number of multiple steps to read in
    if (it + NT_DUMP_NOISE_BUFFER <= it_end) then
      ! full buffer
      nstep_subset_noise_buffer =  NT_DUMP_NOISE_BUFFER
    else
      ! only remaining steps
      nstep_subset_noise_buffer = it_end - it + 1
    endif
    ! checks
    if (nstep_subset_noise_buffer < 1 .or. nstep_subset_noise_buffer > NT_DUMP_NOISE_BUFFER) &
      call exit_MPI(myrank,'Error invalid noise buffer steps')
  endif

  ! index for buffer writing
  ! positive direction (counting forward)
  it_buffer = 1 + icounter_noise_buffer

  ! debug
  !if (myrank == 0) &
  !  print *,'buffer step:',it,'index:',it,'counter:',icounter_noise_buffer,'buffer index:',it_buffer

  ! save current step to buffer
  noise_buffer(:,:,:,:,it_buffer) = noise_surface_movie(:,:,:,:)

  ! writes out buffer (at the end of the cycle)
  if (icounter_noise_buffer == nstep_subset_noise_buffer - 1) then
    ! sets file head index
    ! (incremental forward writing, starting from value 1)
    it_write_index = it - icounter_noise_buffer

    ! writes out buffer to file
    ! debug
    !if (myrank == 0) print *,'  write i/o steps:',nstep_subset_noise_buffer,'index:',it_write_index

    ! reads in buffer
    call write_abs_buffer(9,noise_buffer,reclen_noise * nstep_subset_noise_buffer,it_write_index,reclen_noise)
  endif

  ! updates counter
  if (icounter_noise_buffer < nstep_subset_noise_buffer - 1) then
    icounter_noise_buffer = icounter_noise_buffer + 1
  else
    ! resets counter
    icounter_noise_buffer = 0
  endif

  end subroutine noise_save_surface_movie

!
!-------------------------------------------------------------------------------------------------
!

  subroutine noise_read_add_surface_movie(NGLOB_AB,accel,it_index)

! step 2/3: calculate/reconstruct the "ensemble forward wavefield"
! read surface movie (displacement) at every time steps, injected as the source of "ensemble forward wavefield"
! in step 2, call noise_read_add_surface_movie(..., NSTEP-it+1 ,...)
! in step 3, call noise_read_add_surface_movie(..., it ,...)
!
! there are two subroutines --- noise_read_add_surface_movie_original & noise_read_add_surface_movie
!    noise_read_add_surface_movie_original is implemented at first, which creates one file at each time step
!    noise_read_add_surface_movie is implemented later, which utilizes 'src/shared/binary_c_io.c' for faster I/O,
!                                                   which creates one file for the all time steps
!
! by this modification, the efficiency is greatly improved
! and now, it should be OK to run NOISE_TOMOGRAPHY on a cluster with global storage
!
! called for NOISE_TOMOGRAPHY == 2 and 3

  use specfem_par, only: CUSTOM_REAL,NDIM,NOISE_TOMOGRAPHY,GPU_MODE,UNDO_ATTENUATION, &
    myrank,Mesh_pointer,wgllwgll_xy,it,it_end, &
    NSTEP,NSUBSET_ITERATIONS,NT_DUMP_ATTENUATION,iteration_on_subset,it_subset_end,it_of_this_subset

  use specfem_par_crustmantle, only: NSPEC_TOP,ibelm_top_crust_mantle,ibool_crust_mantle,jacobian2D_top_crust_mantle
  use specfem_par_noise

  implicit none

  integer :: NGLOB_AB
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: accel

  integer,intent(in) :: it_index

  ! local parameters
  integer :: ipoin,ispec2D,ispec,i,j,k,iglob
  real(kind=CUSTOM_REAL) :: eta,fac,jacobianw
  real(kind=CUSTOM_REAL) :: normal_x,normal_y,normal_z

  ! i/o buffer
  integer :: it_read_index,it_buffer
  integer :: it_tmp,it_index_tmp

  ! iteration step
  it_index_tmp = it_index

  ! considers special case for undoing attenuation
  if (UNDO_ATTENUATION) then
    if (NOISE_TOMOGRAPHY == 3) then
      ! this is called from compute_forces_viscoelastic_backward() routine by
      ! > call noise_read_add_surface_movie(NGLOB_CRUST_MANTLE_ADJOINT,b_accel_crust_mantle,it)
      !
      ! time step here should be (reversal of reverse): NSTEP - (NSTEP - it + 1) + 1 = it
      ! for "normal" case of it = 1,..,NSTEP
      !
      ! special case reconstructs wavefield in chunks, going forward from last saved snapshots
      ! example: NT_DUMP_ATTENUATION = 301, NSTEP = 900, NSUBSET_ITERATIONS = 3, iteration_on_subset = 1 -> 3
      !          -> snapshots are at: time it = 0, it = 301 and it = 602
      !          reconstruction goes:
      !              1. subset: it_tmp = 298 down to 1   -> (NSTEP - it_tmp) goes from (602 - 899)
      !              2. subset: it_tmp = 599 down to 299 -> (NSTEP - it_tmp) goes from (301 - 601)
      !              3. subset: it_tmp = 900 down to 600 -> (NSTEP - it_tmp) goes from (0 - 300)
      it_tmp = NSTEP - (NSUBSET_ITERATIONS - iteration_on_subset)*NT_DUMP_ATTENUATION - it_of_this_subset + 1
      ! reversal of reverse:
      !              1. subset: it_index_tmp should go from 298 down to 1
      !              2. subset: it_index_tmp goes from 301 to 601 .. and so on like it_tmp
      it_index_tmp = it_tmp
    endif
  endif

  ! read surface movie
  !
  ! original
  !call read_abs(9,noise_surface_movie,CUSTOM_REAL*NDIM*NGLLX*NGLLY*NSPEC_TOP,it_index)
  !
  ! note: due to heavy file i/o, the reading of the noise surface movie at every timestep becomes very slow
  !       we will be using chunks of multiple time steps
  !
  ! sets buffer steps and read index for file access
  ! (icounter_noise_buffer is a counter between 0 and NT_DUMP_NOISE_BUFFER-1)
  if (icounter_noise_buffer == 0) then
    ! sets number of multiple steps to read in
    if (it + NT_DUMP_NOISE_BUFFER <= it_end) then
      ! full buffer
      nstep_subset_noise_buffer =  NT_DUMP_NOISE_BUFFER
    else
      ! only remaining steps
      nstep_subset_noise_buffer = it_end - it + 1
    endif

    ! special case
    if (NOISE_TOMOGRAPHY == 3 .and. UNDO_ATTENUATION) then
      ! note that buffer size NT_DUMP_NOISE_BUFFER will match attenuation buffer size NT_DUMP_ATTENUATION
      nstep_subset_noise_buffer = it_subset_end
    endif

    ! checks
    if (nstep_subset_noise_buffer < 1 .or. nstep_subset_noise_buffer > NT_DUMP_NOISE_BUFFER) &
      call exit_MPI(myrank,'Error invalid noise buffer steps')
  endif

  ! reads in buffer (at the start of the cycle)
  if (icounter_noise_buffer == 0) then
    ! sets head of index reading
    if (NOISE_TOMOGRAPHY == 2 .or. UNDO_ATTENUATION) then
      ! backward reading, sets head in front to have current it_index_tmp at the end of the buffer array
      it_read_index = it_index_tmp - nstep_subset_noise_buffer + 1
    else if (NOISE_TOMOGRAPHY == 3) then
      ! forward reading
      it_read_index = it_index_tmp
    else
      call exit_MPI(myrank,'Error wrong simulation type for reading noise surface movie')
    endif

    ! debug
    !if (myrank == 0) print *,'  read i/o steps:',nstep_subset_noise_buffer,'index:',it_read_index,it_index_tmp

    ! reads in buffer
    call read_abs_buffer(9,noise_buffer,reclen_noise * nstep_subset_noise_buffer,it_read_index,reclen_noise)
  endif

  ! index for buffer reading
  if (NOISE_TOMOGRAPHY == 2 .or. UNDO_ATTENUATION) then
    ! negative direction (counting back)
    it_buffer = nstep_subset_noise_buffer - icounter_noise_buffer
  else if (NOISE_TOMOGRAPHY == 3) then
    ! positive direction (counting forward)
    it_buffer = 1 + icounter_noise_buffer
  endif

  ! debug
  !if (myrank == 0) &
  !  print *,'buffer step:',it,'index:',it_index_tmp,'counter:',icounter_noise_buffer,'buffer index:',it_buffer

  ! takes current step from buffer
  noise_surface_movie(:,:,:,:) = noise_buffer(:,:,:,:,it_buffer)

  ! updates counter
  if (icounter_noise_buffer < nstep_subset_noise_buffer - 1) then
    icounter_noise_buffer = icounter_noise_buffer + 1
  else
    ! resets counter
    icounter_noise_buffer = 0
  endif

  ! get coordinates of surface mesh and surface displacement
  if (.not. GPU_MODE) then
    ! on CPU
    ipoin = 0
    do ispec2D = 1, NSPEC_TOP ! NSPEC2D_TOP(IREGION_CRUST_MANTLE)
      ispec = ibelm_top_crust_mantle(ispec2D)

      k = NGLLZ

      ! loop on all the points inside the element
      do j = 1,NGLLY
        do i = 1,NGLLX
          ipoin = ipoin + 1
          iglob = ibool_crust_mantle(i,j,k,ispec)

          normal_x = normal_x_noise(ipoin)
          normal_y = normal_y_noise(ipoin)
          normal_z = normal_z_noise(ipoin)

          eta = noise_surface_movie(1,i,j,ispec2D) * normal_x + &
                noise_surface_movie(2,i,j,ispec2D) * normal_y + &
                noise_surface_movie(3,i,j,ispec2D) * normal_z

          ! weighted jacobian
          jacobianw = wgllwgll_xy(i,j) * jacobian2D_top_crust_mantle(i,j,ispec2D)

          ! masks out the noise contributions (e.g. for non-uniform noise distributions)
          fac = eta * jacobianw * mask_noise(ipoin)

          accel(1,iglob) = accel(1,iglob) + fac * normal_x
          accel(2,iglob) = accel(2,iglob) + fac * normal_y
          accel(3,iglob) = accel(3,iglob) + fac * normal_z
        enddo
      enddo

    enddo

  else
    ! on GPU
    call noise_add_surface_movie_gpu(Mesh_pointer,noise_surface_movie)
  endif

  end subroutine noise_read_add_surface_movie

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernels_strength_noise()

! step 3: constructing noise source strength kernel
!
! there are two subroutines --- compute_kernels_strength_noise_original & compute_kernels_strength_noise
!    compute_kernels_strength_noise_original is implemented at first, which creates one file at each time step
!    compute_kernels_strength_noise is implemented later, which utilizes 'src/shared/binary_c_io.c' for faster I/O,
!                                                         which creates only one file for the all time steps
!
! by this modification, the efficiency is greatly improved
! and now, it should be OK to run NOISE_TOMOGRAPHY on a cluster with global storage

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_noise

  implicit none

  ! local parameters
  integer :: i,j,k,ispec,iglob,ipoin,ispec2D
  real(kind=CUSTOM_REAL) :: eta
  real(kind=CUSTOM_REAL) :: normal_x,normal_y,normal_z

  ! read surface movie, needed for sigma_kl_crust_mantle
  !
  ! original
  !call read_abs(9,noise_surface_movie,reclen_noise,it)
  !
  ! note: we read this in for backward call of noise_read_add_surface_movie(),
  !       thus array noise_surface_movie will be filled with values at time it (or retrieved from undo_attenuation buffer)

  ! noise source strength kernel
  ! to keep similar structure to other kernels, the source strength kernel is saved as a volumetric kernel
  ! but only updated at the surface, because the noise is generated there
  if (.not. GPU_MODE) then
    ! on CPU
    ipoin = 0
    do ispec2D = 1, NSPEC_TOP
      ispec = ibelm_top_crust_mantle(ispec2D)

      k = NGLLZ

      ! loop on all the points inside the element
      do j = 1,NGLLY
        do i = 1,NGLLX
          ipoin = ipoin + 1
          iglob = ibool_crust_mantle(i,j,k,ispec)

          normal_x = normal_x_noise(ipoin)
          normal_y = normal_y_noise(ipoin)
          normal_z = normal_z_noise(ipoin)

          eta = noise_surface_movie(1,i,j,ispec2D) * normal_x + &
                noise_surface_movie(2,i,j,ispec2D) * normal_y + &
                noise_surface_movie(3,i,j,ispec2D) * normal_z

          ! daniel: note we might miss the mask here depending on what we want to do;
          !         see Tromp et al. 2010, section 4 and 5:
          !           displ_crust_mantle - the wavefield produced by the 'ensemble adjoint source'
          !           eta                - product \eta_i = G_{jk} \nu_k^{alpha} S_{ij} (eq. 31)
          !         the source region for S_{ij} could be masked for non-uniform noise distributions
          sigma_kl_crust_mantle(i,j,k,ispec) =  sigma_kl_crust_mantle(i,j,k,ispec) &
             + deltat * eta * ( normal_x * displ_crust_mantle(1,iglob) &
                              + normal_y * displ_crust_mantle(2,iglob) &
                              + normal_z * displ_crust_mantle(3,iglob) )
        enddo
      enddo
    enddo

  else
    ! on GPU
    call compute_kernels_strength_noise_gpu(Mesh_pointer,noise_surface_movie,deltat)
  endif

  end subroutine compute_kernels_strength_noise

!
!-------------------------------------------------------------------------------------------------
!

  subroutine save_kernels_strength_noise()

! step 3: save noise source strength kernel

  use specfem_par
  use specfem_par_noise

  implicit none

  ! local parameters
  integer :: ier
  ! local parameters
  real(kind=CUSTOM_REAL) :: scale_kl

  ! scaling factor for kernel units [ s / km^3 ]
  scale_kl = scale_t * scale_displ_inv * 1.d9

  sigma_kl_crust_mantle(:,:,:,:) = sigma_kl_crust_mantle(:,:,:,:) * scale_kl

  ! kernel file output
  if (ADIOS_FOR_KERNELS) then
    call write_kernels_strength_noise_adios()
  else
    ! binary file output
    call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_TMP_PATH)

    open(unit=IOUT_NOISE,file=trim(prname)//'sigma_kernel.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening file sigma_kernel.bin')

    write(IOUT_NOISE) sigma_kl_crust_mantle     ! need to put dimensions back (not done yet)
    close(IOUT_NOISE)
  endif

  end subroutine save_kernels_strength_noise
