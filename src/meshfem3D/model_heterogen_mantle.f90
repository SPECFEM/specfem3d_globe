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
! HMM
!
! generic heterogeneous mantle model
!--------------------------------------------------------------------------------------------------

  module model_heterogen_mantle_par

  use constants, only: CUSTOM_REAL
  implicit none

!
! NOTE: CURRENTLY THIS ROUTINE ONLY WORKS FOR N_R = N_THETA = N_PHI !!!!!
!

  ! heterogen_mantle_model_constants
  integer, parameter :: N_R = 256,N_THETA = 256,N_PHI = 256

  ! model array
  double precision,dimension(:),allocatable :: HMM_rho_in

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: dvpstore

  end module model_heterogen_mantle_par

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_heterogen_mntl_broadcast()

! standard routine to setup model

  use constants, only: myrank,IMAIN,CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,IREGION_CRUST_MANTLE
  use shared_parameters, only: NSPEC_REGIONS

  use model_heterogen_mantle_par

  implicit none

  ! local parameters
  integer :: ier
  integer :: nspec

  ! user info
  if (myrank == 0) then
    write(IMAIN,*) 'broadcast model: heterogen mantle'
    call flush_IMAIN()
  endif

  ! allocates model array
  allocate(HMM_rho_in(N_R*N_THETA*N_PHI),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating HMM array')
  HMM_rho_in(:) = 0._CUSTOM_REAL

  nspec = NSPEC_REGIONS(IREGION_CRUST_MANTLE)
  allocate(dvpstore(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating dvpstore array')
  dvpstore(:,:,:,:) = 0._CUSTOM_REAL

  ! main process reads in model
  if (myrank == 0) then
     write(IMAIN,*) 'Reading in model_heterogen_mantle.'
     call flush_IMAIN()

     call read_heterogen_mantle_model()

     write(IMAIN,*) 'model_heterogen_mantle is read in.'
     call flush_IMAIN()
  endif

  ! broadcast the information read on the main node to all the nodes
  call bcast_all_dp(HMM_rho_in,N_R*N_THETA*N_PHI)

  if (myrank == 0) then
    write(IMAIN,*) 'model_heterogen_mantle is broadcast:'
    write(IMAIN,*) '  First value in HMM:',HMM_rho_in(1)
    write(IMAIN,*) '  Last value in HMM :',HMM_rho_in(N_R*N_THETA*N_PHI)
    call flush_IMAIN()
  endif

  end subroutine model_heterogen_mntl_broadcast

!
!-------------------------------------------------------------------------------------------------
!


  subroutine read_heterogen_mantle_model()

  use constants
  use model_heterogen_mantle_par

  implicit none

  ! local parameters
  integer :: i,num_points,ier

  ! open heterogen.dat
  open(unit=IIN,file='./DATA/heterogen/heterogen.dat',access='direct', &
       form='formatted',recl=20,status='old',action='read',iostat=ier)
  if (ier /= 0 ) call exit_MPI(0,'Error opening model file heterogen.dat')

  ! note: grid is a cubic box around the earth, stretching from [-R_PLANET,+R_PLANET] dimensions
  !
  ! N_R: number of layers in z-direction
  ! N_THETA/N_PHI: horizontal x/y-direction gridding
  num_points = N_R * N_THETA * N_PHI

  do i = 1,num_points
    ! note: single reclength must be 20
    read(IIN,rec=i,fmt='(F20.15)') HMM_rho_in(i)
  enddo

  close(IIN)

  end subroutine read_heterogen_mantle_model

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_heterogen_mantle(ispec,i,j,k,radius,theta,phi,dvs,dvp,drho)

  use constants
  use shared_parameters, only: R_PLANET,THREE_D_MODEL

  use model_heterogen_mantle_par

  implicit none

  ! variable declaration
  integer,intent(in) :: ispec,i,j,k
  double precision,intent(in) :: radius,theta,phi            ! input coordinates
  double precision,intent(out) :: drho,dvp,dvs                ! output anomaly values

  ! local parameters
  double precision :: x,y,z                       ! input converted to Cartesian
  double precision :: x_low,x_high                ! x values used to interpolate
  double precision :: y_low,y_high                ! y values used to interpolate
  double precision :: z_low,z_high                ! z values used to interpolate
  double precision :: delta,delta2                ! weights in record# and in interpolation
  double precision :: rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8 ! rho values at the interpolation points
  double precision :: r_inner,r_outer             ! lower and upper domain bounds for r
  integer :: rec_read                             ! nr of record to be read from heterogen.dat (direct access file)
  double precision :: a,b,c                       ! substitutions in interpolation algorithm (weights)
  double precision :: r_target
  double precision :: lat,lon

  ! modification for heterogeneous model with PREM background reference
  if (THREE_D_MODEL == THREE_D_MODEL_HETEROGEN_PREM) then
    ! chris modif checkers 02/20/21
    lat = (PI/2.0d0-theta)*180.0d0/PI
    lon = phi*180.0d0/PI
    if (lon > 180.0d0) lon = lon - 360.0d0

    ! not fully done yet...
    stop 'MODEL heterogen_prem is not fully implemented yet'

    !call model_heterogen_mantle_prem(ispec,i,j,k,radius,lat,lon,dvs,dvp,drho)

    ! all done
    !return
  endif

  ! dimensions in m
  r_target = radius*R_PLANET

  ! bounds
  ! NOTE: r_outer NEEDS TO BE (just) SMALLER THAN R_PLANET!!!!!!!!
  ! upper/lower radius (in m)
  r_inner = 3.500d6          ! lower bound for heterogeneity zone
  r_outer = R_PLANET - 1.0d1  ! or 6.300d6  !upper bound for heterogeneity zone (lower mantle: e.g. 4.500d6)

  ! increments
  delta = 2.0 * R_PLANET/(real(N_R-1))
  delta2 = 2.0 * R_PLANET/(real(N_R-2))
  !delta2 = 2.0 * R_PLANET/(real(N_R))

  if ((r_target >= r_inner) .and. (r_target <= r_outer)) then
    ! convert spherical point to Cartesian point, move origin to corner
    x = R_PLANET + r_target*sin(theta)*cos(phi)
    y = R_PLANET + r_target*sin(theta)*sin(phi)
    z = R_PLANET + r_target*cos(theta)

    ! determine which points to search for in heterogen.dat
    ! find x_low,y_low,z_low etc.
    x_low = floor(x/delta2) + 1
    x_high = x_low + 1
    y_low = floor(y/delta2) + 1
    y_high = y_low + 1
    z_low = floor(z/delta2) + 1
    z_high = z_low + 1

    ! rho1 at: x_low y_low z_low
    rec_read = int(1+(x_low*N_R*N_R)+(y_low*N_R)+z_low)
    rho1 = HMM_rho_in(rec_read)

    ! rho2 at: x_low y_high z_low
    rec_read = int(1+(x_low*N_R*N_R)+(y_high*N_R)+z_low)
    rho2 = HMM_rho_in(rec_read)

    ! rho3 at: x_high y_low z_low
    rec_read = int(1+(x_high*N_R*N_R)+(y_low*N_R)+z_low)
    rho3 = HMM_rho_in(rec_read)

    ! rho4 at: x_high y_high z_low
    rec_read = int(1+(x_high*N_R*N_R)+(y_high*N_R)+z_low)
    rho4 = HMM_rho_in(rec_read)

    ! rho5 at: x_low y_low z_high
    rec_read = int(1+(x_low*N_R*N_R)+(y_low*N_R)+z_high)
    rho5 = HMM_rho_in(rec_read)

    ! rho6 at: x_low y_high z_high
    rec_read = int(1+(x_low*N_R*N_R)+(y_high*N_R)+z_high)
    rho6 = HMM_rho_in(rec_read)

    ! rho7 at: x_high y_low z_high
    rec_read = int(1+(x_high*N_R*N_R)+(y_low*N_R)+z_high)
    rho7 = HMM_rho_in(rec_read)

    ! rho8 at: x_high y_high z_high
    rec_read = int(1+(x_high*N_R*N_R)+(y_high*N_R)+z_high)
    rho8 = HMM_rho_in(rec_read)

    ! perform linear interpolation between the 8 points
    a = (x-x_low*delta)/delta       ! weight for x
    b = (y-y_low*delta)/delta       ! weight for y
    c = (z-z_low*delta)/delta       ! weight for z

    drho = rho1*(1.-a)*(1.-b)*(1.-c) + rho2*(1.-a)*b*(1.-c) + &
           rho3*a*(1.-b)*(1.-c) + rho4*a*b*(1.-c) + rho5*(1.-a)*(1.-b)*c + &
           rho6*(1.-a)*b*c + rho7*a*(1.-b)*c + rho8*a*b*c

    ! calculate delta vp,vs from the interpolated delta rho
    dvp = (0.55/0.30)*drho
    dvs = (1.00/0.30)*drho

  else
    !outside of heterogeneity domain
    drho = 0.d0
    dvp = 0.d0
    dvs = 0.d0
  endif

  ! stores dvp values for visualizing
  dvpstore(i,j,k,ispec) = real(dvp, kind=CUSTOM_REAL)

  end subroutine model_heterogen_mantle


!
!-------------------------------------------------------------------------------------------------
!



  subroutine model_heterogen_mantle_output_dvp(prname)

  use constants, only: IOUT,MAX_STRING_LEN,myrank
  use model_heterogen_mantle_par

  implicit none

  character(len=MAX_STRING_LEN),intent(in) :: prname

  ! local parameters
  integer :: ier

  ! outputs dvp perturbations
  ! to be visualized with ./xcombine_vol_data util
  open(unit=IOUT,file=prname(1:len_trim(prname))//'dvp.bin', &
        status='unknown',form='unformatted',action='write',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening dvp.bin file')
  write(IOUT) dvpstore
  close(IOUT)

  end subroutine model_heterogen_mantle_output_dvp

! obsolete... only for visualizing needs, which can be done with binary output from above
!
!  ! adios version
!  subroutine model_heterogen_mantle_output_dvp_adios(prname)
!
!  use constants, only: IOUT,MAX_STRING_LEN,myrank
!  use model_heterogen_mantle_par
!  use manager_adios
!  implicit none
!
!  character(len=MAX_STRING_LEN),intent(in) :: prname
!
!  ! local parameters
!  integer :: ier
!
!  write(group_name,"('SPECFEM3D_GLOBE_DVP_reg',i1)") iregion_code
!
!  ! set the adios group size to 0 before incremented by calls to helpers functions.
!  group_size_inc = 0
!  call init_adios_group(myadios_group,group_name)
!
!  !--- Define ADIOS variables -----------------------------
!  local_dim = size (dvpstore)
!  call define_adios_global_array1D(myadios_group, group_size_inc, &
!                                   local_dim, region_name, &
!                                   "dvp", dvpstore)
!  !--- Open an ADIOS handler to the restart file. ---------
!  outputname = trim(LOCAL_PATH) // "/dvp.bp"
!  ! user output
!  if (myrank == 0) write(IMAIN,*) '    saving arrays in ADIOS file: ',trim(outputname)
!
!  if (num_regions_written == 0) then
!    ! opens file for writing
!    call open_file_adios_write(myadios_file,myadios_group,outputname,group_name)
!  else
!    ! opens file for writing in append mode
!    call open_file_adios_write_append(myadios_file,myadios_group,outputname,group_name)
!  endif
!  call set_adios_group_size(myadios_file,group_size_inc)
!
!  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, &
!                                   local_dim, trim(region_name) // "dvp", dvpstore)
!
!  !--- Reset the path to zero and perform the actual write to disk
!  call write_adios_perform(myadios_file)
!  ! closes file
!  call close_file_adios(myadios_file)
!
!  end subroutine model_heterogen_mantle_output_dvp_adios

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_heterogen_mantle_permute_dvp(temp_array_real,perm,nspec)

! permutes dvp array when mesh coloring is desired

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ
  use model_heterogen_mantle_par

  implicit none

  integer,intent(in) :: nspec
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(inout) :: temp_array_real
  integer, dimension(nspec),intent(in) :: perm

  ! local parameters
  integer :: old_ispec,new_ispec

  ! original taken from routine permute_elements_real() in get_perm_color.f90
  ! this copy avoids linking the file for xwrite_profile tool.

  ! copy the original array
  temp_array_real(:,:,:,:) = dvpstore(:,:,:,:)

  do old_ispec = 1,nspec
    new_ispec = perm(old_ispec)
    dvpstore(:,:,:,new_ispec) = temp_array_real(:,:,:,old_ispec)
  enddo

  end subroutine model_heterogen_mantle_permute_dvp
