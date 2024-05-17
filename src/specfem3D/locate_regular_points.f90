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


  subroutine read_kl_regular_grid(GRID)

  use constants, only: IIN,IMAIN,NM_KL_REG_LAYER,PATHNAME_KL_REG, &
    KL_REG_MIN_LAT,KL_REG_MAX_LAT,KL_REG_MIN_LON,KL_REG_MAX_LON

  use specfem_par, only: myrank
  use specfem_par_crustmantle, only: kl_reg_grid_variables
  implicit none

  type (kl_reg_grid_variables), intent(inout) :: GRID

  ! local parameters
  integer :: ios,nlayer,i,nlat,nlon,npts_this_layer
  character(len=256) :: line
  real :: r

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  regular grid input file  : ',trim(PATHNAME_KL_REG)
    call flush_IMAIN()
  endif

  ! improvements to make: read-in by main and broadcast to all secondary processes
  open(IIN,file=trim(PATHNAME_KL_REG),status='old',action='read',iostat=ios)
  if (ios /= 0 ) call exit_MPI(myrank,'Error opening file '//trim(PATHNAME_KL_REG)//' in read_kl_regular_grid() routine')

  ! grid increments
  read(IIN,'(a256)') line
  ! skip comment lines
  do while (line(1:1) == '#')
    read(IIN,'(a256)') line
  enddo

  read(line,*) GRID%dlat, GRID%dlon

  nlayer = 0
  do
    read(IIN,'(a256)',iostat=ios) line
    if (ios /= 0) exit

    ! skip comment lines
    do while (line(1:1) == '#')
      read(IIN,'(a256)',iostat=ios) line
      if (ios /= 0) exit
    enddo

    read(line,*,iostat=ios) r, i
    if (ios /= 0) exit

    if (nlayer >= NM_KL_REG_LAYER) then
      call exit_MPI(myrank, 'Increase NM_KL_REG_LAYER limit')
    endif

    nlayer = nlayer + 1

    GRID%rlayer(nlayer) = r
    GRID%ndoubling(nlayer) = i
  enddo
  close(IIN)

  GRID%nlayer = nlayer

  GRID%npts_total = 0
  GRID%npts_before_layer = 0
  do i = 1, nlayer
    nlon = floor((KL_REG_MAX_LON-KL_REG_MIN_LON)/(GRID%dlon*GRID%ndoubling(i)))+1
    GRID%nlon(i) = nlon
    nlat = floor((KL_REG_MAX_LAT-KL_REG_MIN_LAT)/(GRID%dlat*GRID%ndoubling(i)))+1
    GRID%nlat(i) = nlat
    npts_this_layer = nlon * nlat
    GRID%npts_total = GRID%npts_total + npts_this_layer
    GRID%npts_before_layer(i+1) = GRID%npts_before_layer(i) + npts_this_layer
  enddo
  if (GRID%npts_total <= 0) then
    call exit_MPI(myrank, 'No Model points read in')
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  grid increments (lat/lon): ',GRID%dlat,GRID%dlon
    write(IMAIN,*) '  number of layers         : ',GRID%nlayer
    write(IMAIN,*) '  total number of points   : ',GRID%npts_total
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine read_kl_regular_grid

!==============================================================

  subroutine find_regular_grid_slice_number(slice_number, GRID)

  use constants, only: CUSTOM_REAL,PI,DEGREES_TO_RADIANS, &
    KL_REG_MIN_LAT,KL_REG_MIN_LON

  use specfem_par, only: myrank, addressing, &
    NCHUNKS_VAL, NPROC_XI_VAL, NPROC_ETA_VAL

  use specfem_par_crustmantle, only: kl_reg_grid_variables

  implicit none

  integer, intent(out) :: slice_number(*)

  type (kl_reg_grid_variables), intent(in) :: GRID

  real(kind=CUSTOM_REAL) :: xi_width, eta_width
  integer :: nproc, ilayer, isp, ilat, ilon, k, chunk_isp
  integer :: iproc_xi, iproc_eta
  real :: lat,lon,th,ph,x,y,z,xik,etak,xi_isp,eta_isp,xi1,eta1

  ! assuming 6 chunks full global simulations right now
  if (NCHUNKS_VAL /= 6 .or. NPROC_XI_VAL /= NPROC_ETA_VAL) then
    call exit_MPI(myrank, 'Only deal with 6 chunks at this moment')
  endif

  xi_width = real(PI/2.d0,kind=CUSTOM_REAL); eta_width = real(PI/2.d0,kind=CUSTOM_REAL)
  nproc = NPROC_XI_VAL
  ilayer = 0

  do isp = 1,GRID%npts_total
    if (isp == GRID%npts_before_layer(ilayer+1)+1) ilayer=ilayer+1
    ilat = (isp - GRID%npts_before_layer(ilayer) - 1) / GRID%nlon(ilayer)
    ilon = (isp - GRID%npts_before_layer(ilayer)) - ilat * GRID%nlon(ilayer)

    ! (lat,lon,radius) for isp point
    lat = KL_REG_MIN_LAT + ilat * GRID%dlat * GRID%ndoubling(ilayer)
    th = real((90 - lat) * DEGREES_TO_RADIANS)
    lon = KL_REG_MIN_LON + (ilon - 1) * GRID%dlon * GRID%ndoubling(ilayer)
    ph = real(lon * DEGREES_TO_RADIANS)
    x = sin(th) * cos(ph); y = sin(th) * sin(ph); z = cos(th)

    ! figure out slice number
    chunk_isp = 1; xi_isp = 0; eta_isp = 0
    do k = 1, NCHUNKS_VAL
      call chunk_map(k, x, y, z, xik, etak)
      if (abs(xik) <= PI/4 .and. abs(etak) <= PI/4) then
        chunk_isp = k;  xi_isp = xik; eta_isp = etak; exit
      endif
    enddo
    xi1 = xi_isp / xi_width * 2; eta1 = eta_isp / eta_width * 2
    iproc_xi = floor((xi1+1)/2 * nproc)
    iproc_eta = floor((eta1+1)/2 * nproc)
    slice_number(isp) = addressing(chunk_isp, iproc_xi, iproc_eta)
  enddo

  end subroutine find_regular_grid_slice_number

!==============================================================

! how about using single precision for the iterations?
  subroutine locate_regular_points(npoints_slice_reg,points_slice_reg,GRID, &
                                   nspec,xstore,ystore,zstore,ibool, &
                                   xigll,yigll,zigll,ispec_reg, &
                                   hxir_reg,hetar_reg,hgammar_reg)

  use constants_solver, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NGNOD,NUM_ITER, &
    DEGREES_TO_RADIANS,HUGEVAL,TWO_PI,R_UNIT_SPHERE, &
    NM_KL_REG_PTS,KL_REG_MIN_LAT,KL_REG_MIN_LON

  use shared_parameters, only: R_PLANET

  use specfem_par, only: myrank, &
    NCHUNKS_VAL,NEX_XI_VAL,NEX_ETA_VAL,ANGULAR_WIDTH_XI_IN_DEGREES_VAL,ANGULAR_WIDTH_ETA_IN_DEGREES_VAL

  use specfem_par_crustmantle, only: kl_reg_grid_variables

  implicit none

  ! declarations of regular grid model
  integer, intent(in) :: npoints_slice_reg
  integer, dimension(NM_KL_REG_PTS), intent(in) :: points_slice_reg
  type (kl_reg_grid_variables), intent(in) :: GRID

  ! simulation geometry
  integer, intent(in) :: nspec
  real(kind=CUSTOM_REAL), dimension(*), intent(in) :: xstore,ystore,zstore
  integer, dimension(NGLLX,NGLLY,NGLLZ,*), intent(in) :: ibool

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX), intent(in) :: xigll
  double precision, dimension(NGLLY), intent(in) :: yigll
  double precision, dimension(NGLLZ), intent(in) :: zigll

  ! output
  integer, dimension(NM_KL_REG_PTS), intent(out) :: ispec_reg
  real(kind=CUSTOM_REAL), dimension(NGLLX,NM_KL_REG_PTS), intent(out) :: hxir_reg
  real(kind=CUSTOM_REAL), dimension(NGLLY,NM_KL_REG_PTS), intent(out) :: hetar_reg
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NM_KL_REG_PTS), intent(out) :: hgammar_reg

  ! GLL number of anchors
  integer, dimension(NGNOD) :: anchor_iax,anchor_iay,anchor_iaz

  integer :: i, j, k, isp, ilayer, ilat, ilon, iglob, ix_in, iy_in, iz_in
  integer :: ispec_in, ispec, iter_loop, ia, ipoint
  double precision :: lat, lon, radius, th, ph, x,y,z
  double precision :: x_target, y_target, z_target
  double precision :: distmin_squared,dist_squared
  double precision :: typical_size_squared,element_size
  double precision :: xi,eta,gamma,dx,dy,dz,dxi,deta,dgamma
  double precision :: xix,xiy,xiz
  double precision :: etax,etay,etaz
  double precision :: gammax,gammay,gammaz
  double precision :: ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD
  logical :: locate_target
  double precision, dimension(NGNOD) :: xelm, yelm, zelm

  double precision, dimension(NGLLX) :: hxir
  double precision, dimension(NGLLY) :: hetar
  double precision, dimension(NGLLZ) :: hgammar

  !---------------------------

  call hex_nodes_anchor_ijk(anchor_iax,anchor_iay,anchor_iaz)

  ! compute typical size of elements at the surface
  ! (normalized)
  if (NCHUNKS_VAL == 6) then
    ! estimation for global meshes (assuming 90-degree chunks)
    element_size = TWO_PI * R_UNIT_SPHERE / (4.d0 * NEX_XI_VAL)
  else
    ! estimation for 1-chunk meshes
    ANGULAR_WIDTH_XI_RAD = ANGULAR_WIDTH_XI_IN_DEGREES_VAL * DEGREES_TO_RADIANS
    ANGULAR_WIDTH_ETA_RAD = ANGULAR_WIDTH_ETA_IN_DEGREES_VAL * DEGREES_TO_RADIANS
    element_size = max( ANGULAR_WIDTH_XI_RAD/NEX_XI_VAL,ANGULAR_WIDTH_ETA_RAD/NEX_ETA_VAL ) * R_UNIT_SPHERE
  endif

  ! use 10 times the distance as a criterion for point detections
  typical_size_squared = (10.d0 * element_size)**2

  do ipoint = 1, npoints_slice_reg
    isp = points_slice_reg(ipoint)
    do ilayer = 1, GRID%nlayer
      if (isp <= GRID%npts_before_layer(ilayer+1)) exit
    enddo

    ilat = (isp - GRID%npts_before_layer(ilayer) - 1) / GRID%nlon(ilayer)
    ilon = (isp - GRID%npts_before_layer(ilayer)) - ilat * GRID%nlon(ilayer)

    ! (lat,lon,radius) for isp point
    lat = KL_REG_MIN_LAT + ilat * GRID%dlat * GRID%ndoubling(ilayer)
    lon = KL_REG_MIN_LON + (ilon - 1) * GRID%dlon * GRID%ndoubling(ilayer)

    ! convert radius to meters and then scale
    radius = GRID%rlayer(ilayer) * 1000.0 / R_PLANET

    ! (x,y,z) for isp point
    th = (90.0 - lat) * DEGREES_TO_RADIANS; ph = lon * DEGREES_TO_RADIANS
    x_target = radius * sin(th) * cos(ph)
    y_target = radius * sin(th) * sin(ph)
    z_target = radius * cos(th)

    ! first exclude elements too far away
    locate_target = .false.
    distmin_squared = HUGEVAL
    ix_in = 1
    iy_in = 1
    iz_in = 1
    ispec_in = 1

    do ispec = 1,nspec
      iglob = ibool(1,1,1,ispec)
      dist_squared = (x_target - xstore(iglob))**2 &
                   + (y_target - ystore(iglob))**2 &
                   + (z_target - zstore(iglob))**2

      !  we compare squared distances instead of distances themselves to significantly speed up calculations
      if (dist_squared > typical_size_squared) cycle

      locate_target = .true.
      ! loop only on points inside the element
      ! exclude edges to ensure this point is not
      ! shared with other elements
      ! can be improved if we have a better algorithm of determining if a point
      ! exists inside a 3x3x3 specfem element ???

      do k = 2, NGLLZ-1
        do j = 2, NGLLY-1
          do i = 2, NGLLX-1
            iglob = ibool(i,j,k,ispec)
            dist_squared = (x_target - xstore(iglob))**2 &
                         + (y_target - ystore(iglob))**2 &
                         + (z_target - zstore(iglob))**2
            if (dist_squared < distmin_squared) then
              ix_in = i
              iy_in = j
              iz_in = k
              ispec_in = ispec
              distmin_squared = dist_squared
            endif
          enddo
        enddo
      enddo
    enddo

    if (.not. locate_target) then
      print *, 'Looking for point:', isp, 'layer',ilayer,'ilat/ilon',ilat, ilon
      print *, '  lat/lon = ',lat, lon,'x/y/z = ',x_target, y_target, z_target,'rank', myrank
      call exit_MPI(myrank, 'Error in point_source() array')
    endif

    xi = xigll(ix_in)
    eta = yigll(iy_in)
    gamma = zigll(iz_in)
    ispec_reg(ipoint) = ispec_in

    ! anchors
    do ia = 1, NGNOD
      iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec_in)
      xelm(ia) = dble(xstore(iglob))
      yelm(ia) = dble(ystore(iglob))
      zelm(ia) = dble(zstore(iglob))
    enddo

    ! iterate to solve the nonlinear system
    do iter_loop = 1,NUM_ITER

      ! recompute Jacobian for the new point
      call recompute_jacobian(xelm,yelm,zelm, xi,eta,gamma, x,y,z, &
                              xix,xiy,xiz, etax,etay,etaz, gammax,gammay,gammaz)

      ! compute distance to target location
      dx = - (x - x_target)
      dy = - (y - y_target)
      dz = - (z - z_target)

      ! compute increments
      dxi  = xix*dx + xiy*dy + xiz*dz
      deta = etax*dx + etay*dy + etaz*dz
      dgamma = gammax*dx + gammay*dy + gammaz*dz

      ! update values
      xi = xi + dxi
      eta = eta + deta
      gamma = gamma + dgamma

      ! impose that we stay in that element
      ! (useful if user gives a source outside the mesh for instance)
      if (xi > 1.d0) xi = 1.d0
      if (xi < -1.d0) xi = -1.d0
      if (eta > 1.d0) eta = 1.d0
      if (eta < -1.d0) eta = -1.d0
      if (gamma > 1.d0) gamma = 1.d0
      if (gamma < -1.d0) gamma = -1.d0

    enddo

    ! store l(xi),l(eta),l(gamma)
    call lagrange_any2(xi, NGLLX, xigll, hxir)
    call lagrange_any2(eta, NGLLY, yigll, hetar)
    call lagrange_any2(gamma, NGLLZ, zigll, hgammar)

    hxir_reg(:,ipoint) = real(hxir(:),kind=CUSTOM_REAL)
    hetar_reg(:,ipoint) = real(hetar(:),kind=CUSTOM_REAL)
    hgammar_reg(:,ipoint) = real(hgammar(:),kind=CUSTOM_REAL)

  enddo ! ipoint

  end subroutine locate_regular_points

!==============================================================


  subroutine lagrange_any2(xi,NGLL,xigll,h)

! subroutine to compute the Lagrange interpolants based upon the GLL points
! and their first derivatives at any point xi in [-1,1]

  implicit none

  double precision, intent(in) :: xi
  integer, intent(in) :: NGLL
  double precision, dimension(NGLL), intent(in) :: xigll
  double precision, dimension(NGLL), intent(out) :: h

  integer :: dgr,i
  double precision :: prod1,prod2

  do dgr = 1,NGLL
     prod1 = 1.0d0
     prod2 = 1.0d0

     do i = 1,NGLL
        if (i /= dgr) then
           prod1 = prod1 * (xi         - xigll(i))
           prod2 = prod2 * (xigll(dgr) - xigll(i))
        endif
     enddo

     h(dgr) = prod1 / prod2
  enddo

  end subroutine lagrange_any2

!==============================================================

  subroutine chunk_map(k,xx,yy,zz,xi,eta)

  ! this program get the xi,eta for (xx,yy,zz)
  ! point under the k'th chunk coordinate
  ! transformation

  use constants

  implicit none

  integer, intent(in) :: k
  real, intent(in) :: xx, yy, zz
  real, intent(out) :: xi, eta

  real :: x, y, z
  real, parameter :: EPS = 1e-6

  x = xx; y = yy; z = zz
  if (0 <= x .and. x < EPS)  x = EPS
  if (-EPS < x .and. x < 0)  x = -EPS
  if (0 <= y .and. y < EPS)  y = EPS
  if (-EPS < y .and. y < 0)  y = -EPS
  if (0 <= z .and. z < EPS)  z = EPS
  if (-EPS < z .and. z < 0)  z = -EPS

  if (k == CHUNK_AB) then
     xi = atan(y/z); eta = atan(-x/z)
     if (z < 0)  xi = 10
  else if (k == CHUNK_AC) then
     xi = atan(-z/y); eta = atan(x/y)
     if (y > 0)  xi = 10
  else if (k == CHUNK_BC) then
     xi = atan(-z/x); eta = atan(-y/x)
     if (x > 0)  xi = 10
  else if (k == CHUNK_AC_ANTIPODE) then
     xi = atan(-z/y); eta = atan(-x/y)
     if (y < 0)  xi = 10
  else if (k == CHUNK_BC_ANTIPODE) then
     xi = atan(z/x); eta = atan(-y/x)
     if (x < 0)  xi = 10
  else if (k == CHUNK_AB_ANTIPODE) then
     xi = atan(y/z); eta = atan(x/z)
     if (z > 0)  xi = 10
  else
     stop 'chunk number k < 6'
  endif

  xi = EPS * nint(xi / EPS)
  eta = EPS * nint(eta / EPS)

  end subroutine chunk_map

