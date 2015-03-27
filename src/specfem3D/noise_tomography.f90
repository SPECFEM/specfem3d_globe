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

! subroutine for NOISE TOMOGRAPHY
! characterize noise statistics
! for a given point (xcoord,ycoord,zcoord), specify the noise direction "normal_x/y/z_noise"
!     and noise distribution "mask_noise"
! USERS need to modify this subroutine for their own noise characteristics

  subroutine noise_distribution_direction(xcoord_in,ycoord_in,zcoord_in, &
                  normal_x_noise_out,normal_y_noise_out,normal_z_noise_out, &
                  mask_noise_out)

  use constants

  implicit none

  ! input parameters
  real(kind=CUSTOM_REAL) :: xcoord_in,ycoord_in,zcoord_in
  ! output parameters
  real(kind=CUSTOM_REAL) :: normal_x_noise_out,normal_y_noise_out,normal_z_noise_out,mask_noise_out
  ! local parameters
  real(kind=CUSTOM_REAL) :: xcoord,ycoord,zcoord


  ! coordinates "x/y/zcoord_in" actually contain r theta phi, therefore convert back to x y z
  call rthetaphi_2_xyz(xcoord,ycoord,zcoord, xcoord_in,ycoord_in,zcoord_in)
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

! subroutine for NOISE TOMOGRAPHY
! read parameters

  subroutine read_parameters_noise()

  use specfem_par
  use specfem_par_crustmantle

  implicit none

  ! local parameters
  integer :: ipoin, ispec2D, ispec, i, j, k, iglob, ier
  real(kind=CUSTOM_REAL) :: normal_x_noise_out,normal_y_noise_out,normal_z_noise_out,mask_noise_out
  character(len=MAX_STRING_LEN) :: filename
  real(kind=CUSTOM_REAL), dimension(num_noise_surface_points) :: &
      val_x,val_y,val_z,val_ux,val_uy,val_uz
  real(kind=CUSTOM_REAL), dimension(num_noise_surface_points,0:NPROCTOT_VAL-1) :: &
      val_x_all,val_y_all,val_z_all,val_ux_all,val_uy_all,val_uz_all

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
    open(unit=IOUT_NOISE,file=trim(OUTPUT_FILES)//'/irec_master_noise', &
          status='unknown',action='write',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening output file irec_master_noise')
    write(IOUT_NOISE,*) 'The master receiver is: (RECEIVER ID)', irec_master_noise
    close(IOUT_NOISE)
  endif

  ! checks master irec
  if (irec_master_noise < 1 .or. irec_master_noise > nrec) then
    call exit_MPI(myrank,'Error noise tomography: irec_master_noise is not in range of given number of receivers')
  endif

  ! compute source arrays for "ensemble forward source", which is source of "ensemble forward wavefield"
  if (myrank == islice_selected_rec(irec_master_noise) .OR. myrank == 0) then ! myrank == 0 is used for output only
    call compute_arrays_source_noise(myrank, &
              xi_receiver(irec_master_noise),eta_receiver(irec_master_noise),gamma_receiver(irec_master_noise), &
              nu(:,:,irec_master_noise),noise_sourcearray, xigll,yigll,zigll,NSTEP)
  endif

  ! noise distribution and noise direction
  ipoin = 0

  ! loops over surface points
  ! puts noise distrubution and direction onto the surface points
  do ispec2D = 1, NSPEC_TOP ! NSPEC2D_TOP(IREGION_CRUST_MANTLE)
    ispec = ibelm_top_crust_mantle(ispec2D)

    k = NGLLZ

    ! loop on all the points inside the element
    do j = 1,NGLLY
      do i = 1,NGLLX
        ipoin = ipoin + 1
        iglob = ibool_crust_mantle(i,j,k,ispec)
        ! this subroutine must be modified by USERS
        call noise_distribution_direction(xstore_crust_mantle(iglob), &
                                          ystore_crust_mantle(iglob),zstore_crust_mantle(iglob), &
                                          normal_x_noise_out,normal_y_noise_out,normal_z_noise_out, &
                                          mask_noise_out)

        normal_x_noise(ipoin) = normal_x_noise_out
        normal_y_noise(ipoin) = normal_y_noise_out
        normal_z_noise(ipoin) = normal_z_noise_out
        mask_noise(ipoin)     = mask_noise_out
      enddo
    enddo

  enddo
  ! checks
  if (ipoin /= num_noise_surface_points) call exit_MPI(myrank,'Error invalid number of surface points for noise distribution')

  !!!BEGIN!!! save mask_noise for check, a file called "mask_noise.bin" is saved in "./OUTPUT_FIELS/"
  ipoin = 0
  do ispec2D = 1, NSPEC_TOP ! NSPEC2D_TOP(IREGION_CRUST_MANTLE)
    ispec = ibelm_top_crust_mantle(ispec2D)
    k = NGLLZ
    ! loop on all the points inside the element
    do j = 1,NGLLY
      do i = 1,NGLLX
        ipoin = ipoin + 1
        iglob = ibool_crust_mantle(i,j,k,ispec)
        val_x(ipoin) = xstore_crust_mantle(iglob)
        val_y(ipoin) = ystore_crust_mantle(iglob)
        val_z(ipoin) = zstore_crust_mantle(iglob)
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

  ! save mask_noise data to disk in home directory
  ! this file can be viewed the same way as surface movie data (xcreate_movie_AVS_DX)
  ! create_movie_AVS_DX.f90 needs to be modified in order to do that,
  ! i.e., instead of showing the normal component, change it to either x, y or z component, or the norm.
  if (myrank == 0) then
    open(unit=IOUT_NOISE,file=trim(OUTPUT_FILES)//'/mask_noise.bin', &
              status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening output file mask_noise')

    write(IOUT_NOISE) val_x_all
    write(IOUT_NOISE) val_y_all
    write(IOUT_NOISE) val_z_all
    write(IOUT_NOISE) val_ux_all
    write(IOUT_NOISE) val_uy_all
    write(IOUT_NOISE) val_uz_all

    close(IOUT_NOISE)
  endif
  !!!END!!! save mask_noise for check, a file called "mask_noise.bin" is saved in "./OUTPUT_FIELS/"

  end subroutine read_parameters_noise

!
!-------------------------------------------------------------------------------------------------
!

! subroutine for NOISE TOMOGRAPHY
! check for consistency of the parameters

  subroutine check_parameters_noise()

  use specfem_par
  use specfem_par_crustmantle

  implicit none

  ! local parameters
  integer :: reclen,ier
  integer(kind=8) :: filesize
  character(len=MAX_STRING_LEN) :: outputname

  ! checks if anything to do
  if (NOISE_TOMOGRAPHY == 0) return

  ! info file output
  if (myrank == 0) then
     open(unit=IOUT_NOISE,file=trim(OUTPUT_FILES)//'/NOISE_SIMULATION', &
          status='unknown',action='write',iostat=ier)
     if (ier /= 0 ) call exit_MPI(myrank,'Error opening output file NOISE_SIMULATION')

     write(IOUT_NOISE,*) '*******************************************************************************'
     write(IOUT_NOISE,*) '*******************************************************************************'
     write(IOUT_NOISE,*) 'WARNING!!!!!!!!!!!!'
     write(IOUT_NOISE,*) 'You are running simulations using NOISE TOMOGRAPHY techniques.'
     write(IOUT_NOISE,*) 'Please make sure you understand the procedures before you have a try.'
     write(IOUT_NOISE,*) 'Displacements everywhere at the free surface are saved every timestep,'
     write(IOUT_NOISE,*) 'so make sure that LOCAL_TMP_PATH in DATA/Par_file is not global.'
     write(IOUT_NOISE,*) 'Otherwise the disk storage may be a serious issue, as is the speed of I/O.'
     write(IOUT_NOISE,*) 'Also make sure that NO earthquakes are included,'
     write(IOUT_NOISE,*) 'i.e., set moment tensor to be ZERO in CMTSOLUTION'
     write(IOUT_NOISE,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(IOUT_NOISE,*) 'If you just want a regular EARTHQUAKE simulation,'
     write(IOUT_NOISE,*) 'set NOISE_TOMOGRAPHY = 0 in DATA/Par_file'
     write(IOUT_NOISE,*) '*******************************************************************************'
     write(IOUT_NOISE,*) '*******************************************************************************'
     close(IOUT_NOISE)
  endif

  ! checks parameters
  if (NUMBER_OF_RUNS /= 1 .OR. NUMBER_OF_THIS_RUN /= 1) &
    call exit_mpi(myrank,'NUMBER_OF_RUNS and NUMBER_OF_THIS_RUN must be 1 for NOISE TOMOGRAPHY! check DATA/Par_file')
  if (ROTATE_SEISMOGRAMS_RT) &
    call exit_mpi(myrank,'Do NOT rotate seismograms in the code, change ROTATE_SEISMOGRAMS_RT in DATA/Par_file')
  if (SAVE_ALL_SEISMOS_IN_ONE_FILE .OR. USE_BINARY_FOR_LARGE_FILE) &
    call exit_mpi(myrank,'Please set SAVE_ALL_SEISMOS_IN_ONE_FILE and USE_BINARY_FOR_LARGE_FILE to be .false.')

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

  ! save/read the surface movie using the same c routine as we do for absorbing boundaries (file ID is 9)
  ! size of single record
  reclen=CUSTOM_REAL*NDIM*NGLLX*NGLLY*NSPEC_TOP

  ! check integer size limit: size of b_reclen_field must fit onto an 4-byte integer
  if (NSPEC_TOP > 2147483646 / (CUSTOM_REAL * NGLLX * NGLLY * NDIM)) then
    print *,'reclen of noise surface_movie needed exceeds integer 4-byte limit: ',reclen
    print *,'  ',CUSTOM_REAL, NDIM, NGLLX * NGLLY, NSPEC_TOP
    print*,'bit size fortran: ',bit_size(NSPEC_TOP)
    call exit_MPI(myrank,"error NSPEC_TOP integer limit")
  endif

  ! total file size
  filesize = reclen
  filesize = filesize*NSTEP

  write(outputname,"('/proc',i6.6,'_surface_movie')") myrank
  if (NOISE_TOMOGRAPHY == 1) call open_file_abs_w(9,trim(LOCAL_TMP_PATH)//trim(outputname), &
                                                  len_trim(trim(LOCAL_TMP_PATH)//trim(outputname)), filesize)
  if (NOISE_TOMOGRAPHY == 2) call open_file_abs_r(9,trim(LOCAL_TMP_PATH)//trim(outputname), &
                                                  len_trim(trim(LOCAL_TMP_PATH)//trim(outputname)), filesize)
  if (NOISE_TOMOGRAPHY == 3) call open_file_abs_r(9,trim(LOCAL_TMP_PATH)//trim(outputname), &
                                                  len_trim(trim(LOCAL_TMP_PATH)//trim(outputname)), filesize)

  end subroutine check_parameters_noise

!
!-------------------------------------------------------------------------------------------------
!

! subroutine for NOISE TOMOGRAPHY
! read and construct the "source" (source time function based upon noise spectrum)
! for "ensemble forward source"

  subroutine compute_arrays_source_noise(myrank, &
                                         xi_noise,eta_noise,gamma_noise,nu_single,noise_sourcearray, &
                                         xigll,yigll,zigll,NSTEP)

  use constants_solver
  use shared_parameters,only:OUTPUT_FILES,DT

  implicit none

  ! input parameters
  integer :: myrank, NSTEP
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
  !double precision,parameter :: scale_displ_inv = 1.d0/R_EARTH ! non-dimensional scaling
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
      print*,'Error noise source S_squared file length: NSTEP length required is ',NSTEP,' with time step size ',DT
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
      print*,'Error noise nu_master file length: number of required components is 3'
      call exit_MPI(myrank,&
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

! subroutine for NOISE TOMOGRAPHY
! step 1: calculate the "ensemble forward source"
! add noise spectrum to the location of master receiver

  subroutine noise_add_source_master_rec()

! the first step of noise tomography is to use |S(\omega)|^2 as a point force source at one of the receivers.
! hence, instead of a moment tensor 'sourcearrays', a 'noise_sourcearray' for a point force is needed.
! furthermore, the CMTSOLUTION needs to be zero, i.e., no earthquakes.
! now this must be manually set in DATA/CMTSOLUTION, by USERS.

  use specfem_par
  use specfem_par_crustmantle

  implicit none

  ! local parameters
  integer :: i,j,k,iglob

  ! adds noise source (only if this proc carries the noise)
  if (.not. GPU_MODE) then
    ! on CPU
    if (myrank == islice_selected_rec(irec_master_noise)) then
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
    endif

  else
    ! on GPU
    call noise_add_source_master_rec_gpu(Mesh_pointer,it,irec_master_noise,islice_selected_rec)
  endif

  end subroutine noise_add_source_master_rec

!
!-------------------------------------------------------------------------------------------------
!

! subroutine for NOISE TOMOGRAPHY
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

  subroutine noise_save_surface_movie()

  use specfem_par
  use specfem_par_crustmantle

  implicit none

  ! local parameters
  integer :: ispec2D,ispec,i,j,k,iglob

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
  call write_abs(9,noise_surface_movie,CUSTOM_REAL*NDIM*NGLLX*NGLLY*NSPEC_TOP,it)

  end subroutine noise_save_surface_movie

!
!-------------------------------------------------------------------------------------------------
!

! subroutine for NOISE TOMOGRAPHY
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

  subroutine noise_read_add_surface_movie(NGLOB_AB,accel,it_index)

  use specfem_par
  use specfem_par_crustmantle

  implicit none

  integer :: NGLOB_AB
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB),intent(inout) :: accel

  integer,intent(in) :: it_index

  ! local parameters
  integer :: ipoin,ispec2D,ispec,i,j,k,iglob
  real(kind=CUSTOM_REAL) :: eta

  ! read surface movie
  call read_abs(9,noise_surface_movie,CUSTOM_REAL*NDIM*NGLLX*NGLLY*NSPEC_TOP,it_index)

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

          eta = noise_surface_movie(1,i,j,ispec2D) * normal_x_noise(ipoin) + &
                noise_surface_movie(2,i,j,ispec2D) * normal_y_noise(ipoin) + &
                noise_surface_movie(3,i,j,ispec2D) * normal_z_noise(ipoin)

          accel(1,iglob) = accel(1,iglob) + eta * mask_noise(ipoin) * normal_x_noise(ipoin) &
                              * wgllwgll_xy(i,j) * jacobian2D_top_crust_mantle(i,j,ispec2D)
          accel(2,iglob) = accel(2,iglob) + eta * mask_noise(ipoin) * normal_y_noise(ipoin) &
                              * wgllwgll_xy(i,j) * jacobian2D_top_crust_mantle(i,j,ispec2D)
          accel(3,iglob) = accel(3,iglob) + eta * mask_noise(ipoin) * normal_z_noise(ipoin) &
                              * wgllwgll_xy(i,j) * jacobian2D_top_crust_mantle(i,j,ispec2D)
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

! subroutine for NOISE TOMOGRAPHY
! step 3: constructing noise source strength kernel
!
! there are two subroutines --- compute_kernels_strength_noise_original & compute_kernels_strength_noise
!    compute_kernels_strength_noise_original is implemented at first, which creates one file at each time step
!    compute_kernels_strength_noise is implemented later, which utilizes 'src/shared/binary_c_io.c' for faster I/O,
!                                                         which creates only one file for the all time steps
!
! by this modification, the efficiency is greatly improved
! and now, it should be OK to run NOISE_TOMOGRAPHY on a cluster with global storage

  subroutine compute_kernels_strength_noise()

  use specfem_par
  use specfem_par_crustmantle

  implicit none

  ! local parameters
  integer :: i,j,k,ispec,iglob,ipoin,ispec2D
  real(kind=CUSTOM_REAL) :: eta

  ! read surface movie, needed for sigma_kl_crust_mantle
  call read_abs(9,noise_surface_movie,CUSTOM_REAL*NDIM*NGLLX*NGLLY*NSPEC_TOP,it)

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

          eta = noise_surface_movie(1,i,j,ispec2D) * normal_x_noise(ipoin) + &
                noise_surface_movie(2,i,j,ispec2D) * normal_y_noise(ipoin) + &
                noise_surface_movie(3,i,j,ispec2D) * normal_z_noise(ipoin)

          sigma_kl_crust_mantle(i,j,k,ispec) =  sigma_kl_crust_mantle(i,j,k,ispec) &
             + deltat * eta * ( normal_x_noise(ipoin) * displ_crust_mantle(1,iglob) &
                              + normal_y_noise(ipoin) * displ_crust_mantle(2,iglob) &
                              + normal_z_noise(ipoin) * displ_crust_mantle(3,iglob) )
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

! subroutine for NOISE TOMOGRAPHY
! step 3: save noise source strength kernel

  subroutine save_kernels_strength_noise()

  use specfem_par
  use specfem_par_crustmantle

  implicit none

  ! local parameters
  integer :: ier
  ! local parameters
  real(kind=CUSTOM_REAL) :: scale_kl

  ! scaling factor for kernel units [ s / km^3 ]
  scale_kl = scale_t/scale_displ * 1.d9

  sigma_kl_crust_mantle(:,:,:,:) = sigma_kl_crust_mantle(:,:,:,:) * scale_kl

  ! kernel file output
  if (ADIOS_FOR_KERNELS) then
    call write_kernels_strength_noise_adios(current_adios_handle)
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
