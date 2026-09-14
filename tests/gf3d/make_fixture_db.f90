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

!----
!---- make_fixture_db -- a synthetic Green function database, ~4.5 MB
!----
!---- Writes a database that gf_open accepts and that gf_locate_source and
!---- the extraction path can run against, so that the database-backed tests
!---- run on a bare checkout instead of needing the 786 MB gitignored example.
!---- It is built into a temporary directory and removed again, never committed.
!---- tests/gf3d/gfdb_env.sh builds one into a temporary directory whenever
!---- GF3D_TEST_GFDB is unset.
!----
!---- usage: make_fixture_db <directory>       writes <directory>/GFDB/...
!----
!---- What this is NOT
!---- ----------------
!---- It is not an oracle. The displacement it writes is a smooth analytic
!---- function chosen to be finite, causal and reproducible, not a Green
!---- function; no test may assert a physical value against it. The solver's
!---- forward run remains the only oracle for anything physical, exactly as
!---- gf3df_integration_plan/testing.md's coverage table says. What the
!---- fixture pins is everything structural: that gf_open reads the layout,
!---- that handles are independent, that the kd-tree is owned correctly, that
!---- the C ABI and the ctypes binding agree, and that the error paths fire.
!----
!---- Why the elements are affine
!---- --------------------------
!---- gf_check_anchors pulls the 27 anchors out of the stored xyz(3,5,5,5) and
!---- re-evaluates the tri-quadratic shape functions at every GLL point,
!---- requiring agreement to GF_ANCHOR_TOL = 1e-6. A tri-quadratic reproduces
!---- an affine map exactly, so an affine hexahedron passes at the float32
!---- storage floor (~1e-8) with three decades to spare. A curved *radial*
!---- element -- each GLL node pushed onto its own sphere -- is NOT
!---- tri-quadratic and trips the guard exactly as a USE_GLL = .true.
!---- database does, so do not be tempted to make this look more like a real
!---- mesh than it is.
!----
!---- The elements are placed as a Cartesian box around the position of
!---- REF_DATA/CMTSOLUTION under a plain spherical conversion. That placement
!---- is deliberately coarse -- the box half-width is ~300 km, against a
!---- geographic-to-geocentric latitude difference of at most ~21 km -- so
!---- that no test depends on this file reproducing the library's geographic
!---- chain. It only has to put the point somewhere inside the element.
!----

  program make_fixture_db

  use hdf5

  implicit none

  ! the build's own values; a database whose ngll differs is refused by
  ! gf_validate, so there is nothing to choose here
  integer, parameter :: NG = 5
  integer, parameter :: NCOMP = 3

  ! nt_subsampled must be exactly nstep/subsample_step.
  !
  ! 2048 * 0.1 s = 205 s of record. It could be far shorter and still be a
  ! valid database, but REF_DATA/CMTSOLUTION has hdur = 60 s and so asks for
  ! t0 = 90 s: on a 6 s axis the output would be 93 % prepended zeros, and a
  ! test comparing traces would be comparing padding. This keeps the stored
  ! part the majority of the trace. It costs 4.6 MB in a temp directory.
  integer, parameter :: NSTEP = 2048
  integer, parameter :: SUBSAMPLE_STEP = 8
  integer, parameter :: NT_SUB = NSTEP / SUBSAMPLE_STEP

  double precision, parameter :: DT = 0.1d0
  double precision, parameter :: T0 = 3.0d0
  double precision, parameter :: R_PLANET = 6371000.d0
  double precision, parameter :: RHOAV = 5514.3d0
  double precision, parameter :: SCALE_DISPL = 1.d0

  ! REF_DATA/CMTSOLUTION's position, and a box wide enough that the exact
  ! latitude convention cannot move the point out of it
  double precision, parameter :: SRC_LAT = -5.8120d0
  double precision, parameter :: SRC_LON = -75.2700d0
  double precision, parameter :: SRC_DEPTH_KM = 122.6d0
  double precision, parameter :: HALF_WIDTH = 0.05d0        ! ~320 km, non-dimensional

  integer, parameter :: NELEM = 2
  integer, parameter :: NSTA = 2

  character(len=8), parameter :: NET(NSTA) = (/ 'II      ', 'IU      ' /)
  character(len=8), parameter :: STA(NSTA) = (/ 'FIX1    ', 'FIX2    ' /)

  character(len=512) :: root,gfdb,edir,fname
  character(len=16) :: hexcode(NELEM)
  double precision :: centre(3),cen(3,NELEM)
  double precision :: xyz(3,NG,NG,NG)
  real :: xyz32(3,NG,NG,NG)
  real :: displ(NCOMP,NCOMP,NG,NG,NG,NT_SUB)
  real :: stf(NSTEP)
  integer :: ie,is,ierr,narg

  narg = command_argument_count()
  if (narg < 1) then
    write(*,'(a)') 'usage: make_fixture_db <directory>'
    stop 1
  endif
  call get_command_argument(1,root)

  gfdb = trim(root)//'/GFDB'

  call h5open_f(ierr)
  if (ierr /= 0) call die('could not initialise HDF5')

  call make_dir(trim(gfdb))
  call make_dir(trim(gfdb)//'/elements')
  call make_dir(trim(gfdb)//'/stations')

  !--- where to put the elements ---------------------------------------

  call sphere_position(SRC_LAT,SRC_LON,SRC_DEPTH_KM,centre)

  ! Two elements side by side along x, so the index has more than one entry
  ! and the Morton codes have something to order. The codes are made up: the
  ! library only requires 16 uppercase hex digits, strictly ascending.
  hexcode(1) = '0000000000000010'
  hexcode(2) = '00000000000000A0'

  !--- mesh_info.h5 -----------------------------------------------------

  call write_mesh_info(trim(gfdb)//'/mesh_info.h5')

  !--- the elements -----------------------------------------------------

  do ie = 1,NELEM
    edir = trim(gfdb)//'/elements/'//trim(hexcode(ie))
    call make_dir(trim(edir))

    call affine_element(centre,dble(ie-1)*HALF_WIDTH,HALF_WIDTH,xyz,cen(:,ie))
    xyz32(:,:,:,:) = real(xyz(:,:,:,:))

    call write_coordinates(trim(edir)//'/coordinates.h5',xyz32,cen(:,ie))

    do is = 1,NSTA
      call fill_displacement(xyz,cen(:,ie),is,displ)
      fname = trim(edir)//'/'//trim(NET(is))//'.'//trim(STA(is))//'.h5'
      call write_element_station(trim(fname),displ)
    enddo
  enddo

  !--- the stations -----------------------------------------------------

  do is = 1,NSTA
    call fill_stf(stf)
    fname = trim(gfdb)//'/stations/'//trim(NET(is))//'.'//trim(STA(is))//'.h5'
    call write_station(trim(fname),is,stf)
  enddo

  call h5close_f(ierr)

  write(*,'(a,a)') 'fixture database written to ',trim(gfdb)
  write(*,'(a,i0,a,i0,a,i0,a)') '  ',NELEM,' elements, ',NSTA,' stations, ',NT_SUB,' samples'

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine die(msg)

  implicit none
  character(len=*), intent(in) :: msg

  write(*,'(a)') 'make_fixture_db: '//msg
  stop 1

  end subroutine die

!
!-------------------------------------------------------------------------------------------------
!

  subroutine make_dir(path)

! mkdir -p. execute_command_line rather than a C shim: this is a test
! generator, it runs once, and the shim would be one more thing to build.

  implicit none
  character(len=*), intent(in) :: path

  integer :: stat

  call execute_command_line('mkdir -p '''//path//'''',exitstat=stat)
  if (stat /= 0) call die('could not create '//path)

  end subroutine make_dir

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sphere_position(lat,lon,depth_km,xyz_out)

! A plain spherical geographic-to-Cartesian conversion, non-dimensionalised
! by R_PLANET. Used ONLY to decide where to put the element -- see the header:
! nothing asserts against it, and it deliberately ignores the geocentric
! latitude correction the library applies.

  implicit none
  double precision, intent(in) :: lat,lon,depth_km
  double precision, dimension(3), intent(out) :: xyz_out

  double precision, parameter :: PI_ = 3.141592653589793d0
  double precision :: theta,phi,r

  theta = (90.d0 - lat)*PI_/180.d0
  phi = lon*PI_/180.d0
  r = 1.d0 - depth_km*1000.d0/R_PLANET

  xyz_out(1) = r*sin(theta)*cos(phi)
  xyz_out(2) = r*sin(theta)*sin(phi)
  xyz_out(3) = r*cos(theta)

  end subroutine sphere_position

!
!-------------------------------------------------------------------------------------------------
!

  subroutine affine_element(centre_in,offset,h,xyz_out,centroid)

! an affine hexahedron: x(i,j,k) = c + A (xi_i, eta_j, gam_k)
!
! A is non-symmetric and its columns are linearly independent, so the
! Jacobian is constant and positive and gf_find_local_coords converges in one
! Newton step. Affine is what makes the 27-anchor guard exact -- see header.

  implicit none
  double precision, dimension(3), intent(in) :: centre_in
  double precision, intent(in) :: offset,h
  double precision, dimension(3,NG,NG,NG), intent(out) :: xyz_out
  double precision, dimension(3), intent(out) :: centroid

  double precision, dimension(NG) :: xi
  double precision, dimension(3,3) :: amat
  double precision, dimension(3) :: c
  integer :: i,j,k,p

  call gll_points(xi)

  c(:) = centre_in(:)
  c(1) = c(1) + offset

  amat(1,1) = h        ; amat(1,2) = 0.12d0*h ; amat(1,3) = -0.05d0*h
  amat(2,1) = -0.08d0*h ; amat(2,2) = h        ; amat(2,3) = 0.10d0*h
  amat(3,1) = 0.06d0*h  ; amat(3,2) = -0.09d0*h ; amat(3,3) = h

  do k = 1,NG
    do j = 1,NG
      do i = 1,NG
        do p = 1,3
          xyz_out(p,i,j,k) = c(p) + amat(p,1)*xi(i) + amat(p,2)*xi(j) + amat(p,3)*xi(k)
        enddo
      enddo
    enddo
  enddo

  centroid(:) = c(:)

  end subroutine affine_element

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gll_points(xi)

! the 5-point Gauss-Lobatto-Legendre abscissae on [-1,1]
!
! Written out rather than obtained from zwgljd: the anchors only have to be
! consistent with the stored coordinates, and an affine map is reproduced
! exactly at any node placement, so this file need not link the solver's GLL
! machinery to be correct.

  implicit none
  double precision, dimension(NG), intent(out) :: xi

  xi(1) = -1.d0
  xi(2) = -sqrt(3.d0/7.d0)
  xi(3) = 0.d0
  xi(4) = sqrt(3.d0/7.d0)
  xi(5) = 1.d0

  end subroutine gll_points

!
!-------------------------------------------------------------------------------------------------
!

  subroutine fill_displacement(xyz_in,centroid,ista,d)

! a smooth, finite, causal analytic field -- not a Green function
!
! Indices are the reader's: d(force component, displacement component, i,j,k,
! time), see gf_element_io.F90's header.
!
! Causal on purpose. gf_seismograms warns when a trace is not silent before
! the reciprocal source switches on, because the conversion zero-extends
! there and would be in error -- a real Green function has nothing before the
! first arrival. A field that is non-zero at t = 0 makes every extraction
! from this fixture print that warning, which would read as a defect in the
! library rather than a property of the fixture. So the field is exactly zero
! until T_ONSET and leaves it with a C2 ramp.

  implicit none
  double precision, dimension(3,NG,NG,NG), intent(in) :: xyz_in
  double precision, dimension(3), intent(in) :: centroid
  integer, intent(in) :: ista
  real, dimension(NCOMP,NCOMP,NG,NG,NG,NT_SUB), intent(out) :: d

  ! a fifth of the record, comfortably past the kernel half-length
  double precision, parameter :: T_ONSET = 0.2d0*dble(NSTEP)*DT

  ! Chosen so the *extracted seismogram* lands around 1e-6 m, the range the
  ! shipped example produces. The contraction divides the file's dyne-cm by
  ! scale_moment ~ 1e28, so a physically plausible displacement here would put
  ! the output at ~1e-36 -- a decade and a half above the float32 subnormal
  ! floor, which is precisely the band where a writer's real(x) under -ftz
  ! flushes to zero and a bitwise SAC comparison fails on one compiler and not
  ! another. The number is meaningless as a displacement; the fixture is not
  ! physical. It is chosen to keep the arithmetic in a range the tests tuned
  ! for the real database already handle.
  double precision, parameter :: AMPLITUDE = 1.d18

  integer :: i,j,k,a,b,it
  double precision :: t,s,env,tau

  do it = 1,NT_SUB
    t = dble(it-1)*DT*dble(SUBSAMPLE_STEP)

    if (t <= T_ONSET) then
      env = 0.d0
    else
      tau = t - T_ONSET
      ! tau**3 makes value, slope and curvature vanish at the onset, so the
      ! trace is smooth where it starts rather than stepping
      env = tau**3 * exp(-0.05d0*tau) * sin(0.35d0*tau + 0.4d0*dble(ista))
    endif

    do k = 1,NG
      do j = 1,NG
        do i = 1,NG
          s = shape_factor(xyz_in(:,i,j,k),centroid)
          do b = 1,NCOMP
            do a = 1,NCOMP
              d(a,b,i,j,k,it) = real(AMPLITUDE*env*(1.d0 + 0.1d0*dble(a) + 0.03d0*dble(b))*s)
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  end subroutine fill_displacement

!
!-------------------------------------------------------------------------------------------------
!

  double precision function shape_factor(x,centroid)

! how the field varies across the element, in element-local coordinates
!
! A *cubic*, and that is the point. test_gf_partials_db checks that a central
! difference of the seismogram against position converges at second order,
! by taking the ratio of the errors at its two largest steps and requiring
! it in [3,5]. The truncation error of a central difference is (h^2/6) f''',
! so a field with no third derivative -- anything linear in position -- has
! no truncation error to measure, the finite difference sees only the float32
! storage noise, and the ratio comes out at 1. A linear field was what this
! fixture had first, and that check failed on it.
!
! u is scaled by the element half-width, so the field varies by O(1) across
! one element and f'''/f' ~ (1/HALF_WIDTH)^2 ~ 1.6e3. At the test's largest
! step (4e-3 degrees, ~7e-5 non-dimensional) that puts the relative
! truncation error near 1e-6, two decades clear of the noise floor, and it
! falls by four for each halving of the step, as the assertion expects.
!
! Degree 3 in each variable is inside the degree-4 GLL basis, so what is
! stored at the nodes is reproduced exactly between them: the interpolant the
! library differentiates is this polynomial and not an approximation of it.

  implicit none
  double precision, dimension(3), intent(in) :: x,centroid

  double precision :: u1,u2,u3

  u1 = (x(1) - centroid(1))/HALF_WIDTH
  u2 = (x(2) - centroid(2))/HALF_WIDTH
  u3 = (x(3) - centroid(3))/HALF_WIDTH

  ! u2 carries a strong cubic and a weak linear term on purpose. Depth moves
  ! the point along the radius, which at this source position is almost
  ! entirely -y, i.e. almost entirely u2; and the convergence check divides by
  ! the partial itself, so a large f''' over a small f' is what lifts the
  ! depth truncation clear of the storage noise. With u2 shaped like the other
  ! two the depth ratio came out at 1.4, not 4.
  shape_factor = 1.d0 &
               + 0.60d0*u1 + 0.50d0*u1*u1 + 0.40d0*u1*u1*u1 &
               + 0.05d0*u2 + 0.35d0*u2*u2 + 1.20d0*u2*u2*u2 &
               + 0.25d0*u3 + 0.22d0*u3*u3 + 0.30d0*u3*u3*u3 &
               + 0.20d0*u1*u2 - 0.15d0*u2*u3

  end function shape_factor

!
!-------------------------------------------------------------------------------------------------
!

  subroutine fill_stf(s)

! the reciprocal run's source time function, on the *unsubsampled* axis:
! gf_read_stations requires exactly nstep samples

  implicit none
  real, dimension(NSTEP), intent(out) :: s

  integer :: it
  double precision :: t,hd

  hd = 4.d0*DT
  do it = 1,NSTEP
    t = dble(it-1)*DT - T0
    s(it) = real(exp(-(t/hd)**2))
  enddo

  end subroutine fill_stf

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_mesh_info(path)

! the root attributes gf_read_mesh_info reads
!
! TOPOGRAPHY and ELLIPTICITY are 0, so there are no datasets at all: no
! ibathy_topo (which is 58 MB of the real example's mesh_info.h5) and no
! ellipticity splines. The flags must still be present -- gf_read_flag errors
! on a missing attribute, it does not default.

  implicit none
  character(len=*), intent(in) :: path

  integer(kind=HID_T) :: fid
  integer :: ier

  call h5fcreate_f(path,H5F_ACC_TRUNC_F,fid,ier)
  if (ier /= 0) call die('could not create '//path)

  call put_attr_d(fid,'dt',DT)
  call put_attr_d(fid,'t0',T0)
  call put_attr_i(fid,'nstep',NSTEP)
  call put_attr_i(fid,'nt_subsampled',NT_SUB)
  call put_attr_i(fid,'subsample_step',SUBSAMPLE_STEP)
  call put_attr_i(fid,'buffer_size',1)
  call put_attr_i(fid,'neighbor_shells',0)
  call put_attr_i(fid,'ngllx',NG)
  call put_attr_i(fid,'nglly',NG)
  call put_attr_i(fid,'ngllz',NG)
  call put_attr_d(fid,'scale_displ',SCALE_DISPL)
  call put_attr_d(fid,'R_PLANET',R_PLANET)
  call put_attr_d(fid,'RHOAV',RHOAV)
  call put_attr_i(fid,'TOPOGRAPHY',0)
  call put_attr_i(fid,'ELLIPTICITY',0)
  call put_attr_i(fid,'rotation',0)
  call put_attr_i(fid,'attenuation',0)
  call put_attr_i(fid,'gravity',0)

  call h5fclose_f(fid,ier)

  end subroutine write_mesh_info

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_coordinates(path,xyz32_in,centroid)

  implicit none
  character(len=*), intent(in) :: path
  real, dimension(3,NG,NG,NG), intent(in) :: xyz32_in
  double precision, dimension(3), intent(in) :: centroid

  integer(kind=HID_T) :: fid,sid,did
  integer(kind=HSIZE_T) :: dims(4)
  integer :: ier

  call h5fcreate_f(path,H5F_ACC_TRUNC_F,fid,ier)
  if (ier /= 0) call die('could not create '//path)

  dims = (/ 3_HSIZE_T, int(NG,HSIZE_T), int(NG,HSIZE_T), int(NG,HSIZE_T) /)
  call h5screate_simple_f(4,dims,sid,ier)
  call h5dcreate_f(fid,'xyz',H5T_NATIVE_REAL,sid,did,ier)
  call h5dwrite_f(did,H5T_NATIVE_REAL,xyz32_in,dims,ier)
  call h5dclose_f(did,ier)
  call h5sclose_f(sid,ier)

  ! read by the directory-scan index route (gf_scan_element_dirs)
  call put_attr_d(fid,'cx',centroid(1))
  call put_attr_d(fid,'cy',centroid(2))
  call put_attr_d(fid,'cz',centroid(3))

  call h5fclose_f(fid,ier)

  end subroutine write_coordinates

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_element_station(path,d)

  implicit none
  character(len=*), intent(in) :: path
  real, dimension(NCOMP,NCOMP,NG,NG,NG,NT_SUB), intent(in) :: d

  integer(kind=HID_T) :: fid,sid,did
  integer(kind=HSIZE_T) :: dims(6)
  integer :: ier

  call h5fcreate_f(path,H5F_ACC_TRUNC_F,fid,ier)
  if (ier /= 0) call die('could not create '//path)

  dims = (/ int(NCOMP,HSIZE_T), int(NCOMP,HSIZE_T), int(NG,HSIZE_T), &
            int(NG,HSIZE_T), int(NG,HSIZE_T), int(NT_SUB,HSIZE_T) /)
  call h5screate_simple_f(6,dims,sid,ier)
  call h5dcreate_f(fid,'displacement',H5T_NATIVE_REAL,sid,did,ier)
  call h5dwrite_f(did,H5T_NATIVE_REAL,d,dims,ier)
  call h5dclose_f(did,ier)
  call h5sclose_f(sid,ier)

  ! the only completion flag gf_check_completion inspects
  call put_attr_i(fid,'computed_ALL',1)

  call h5fclose_f(fid,ier)

  end subroutine write_element_station

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_station(path,ista,s)

  implicit none
  character(len=*), intent(in) :: path
  integer, intent(in) :: ista
  real, dimension(NSTEP), intent(in) :: s

  integer(kind=HID_T) :: fid,sid,did
  integer(kind=HSIZE_T) :: dims(1)
  integer :: ier

  call h5fcreate_f(path,H5F_ACC_TRUNC_F,fid,ier)
  if (ier /= 0) call die('could not create '//path)

  ! all seven are required by gf_read_stations
  call put_attr_d(fid,'latitude',10.d0 + 5.d0*dble(ista))
  call put_attr_d(fid,'longitude',-70.d0 + 3.d0*dble(ista))
  call put_attr_d(fid,'depth',0.d0)
  call put_attr_d(fid,'hdur',4.d0*DT)
  call put_attr_d(fid,'f_cutoff',0.d0)
  call put_attr_d(fid,'factor_force_source',1.d15)
  call put_attr_d(fid,'time_shift',0.d0)

  ! nstep samples, not nt_subsampled: gf_read_stations checks this
  dims = (/ int(NSTEP,HSIZE_T) /)
  call h5screate_simple_f(1,dims,sid,ier)
  call h5dcreate_f(fid,'stf',H5T_NATIVE_REAL,sid,did,ier)
  call h5dwrite_f(did,H5T_NATIVE_REAL,s,dims,ier)
  call h5dclose_f(did,ier)
  call h5sclose_f(sid,ier)

  call h5fclose_f(fid,ier)

  end subroutine write_station

!
!-------------------------------------------------------------------------------------------------
!

  subroutine put_attr_d(loc,name,val)

! a scalar written as a 1-element array, which is the shape the writer uses
! and therefore the shape the reader expects

  implicit none
  integer(kind=HID_T), intent(in) :: loc
  character(len=*), intent(in) :: name
  double precision, intent(in) :: val

  integer(kind=HID_T) :: sid,aid
  integer(kind=HSIZE_T) :: dims(1)
  double precision :: buf(1)
  integer :: ier

  dims = (/ 1_HSIZE_T /)
  buf(1) = val
  call h5screate_simple_f(1,dims,sid,ier)
  call h5acreate_f(loc,name,H5T_NATIVE_DOUBLE,sid,aid,ier)
  call h5awrite_f(aid,H5T_NATIVE_DOUBLE,buf,dims,ier)
  call h5aclose_f(aid,ier)
  call h5sclose_f(sid,ier)

  end subroutine put_attr_d

!
!-------------------------------------------------------------------------------------------------
!

  subroutine put_attr_i(loc,name,val)

  implicit none
  integer(kind=HID_T), intent(in) :: loc
  character(len=*), intent(in) :: name
  integer, intent(in) :: val

  integer(kind=HID_T) :: sid,aid
  integer(kind=HSIZE_T) :: dims(1)
  integer :: buf(1)
  integer :: ier

  dims = (/ 1_HSIZE_T /)
  buf(1) = val
  call h5screate_simple_f(1,dims,sid,ier)
  call h5acreate_f(loc,name,H5T_NATIVE_INTEGER,sid,aid,ier)
  call h5awrite_f(aid,H5T_NATIVE_INTEGER,buf,dims,ier)
  call h5aclose_f(aid,ier)
  call h5sclose_f(sid,ier)

  end subroutine put_attr_i

  end program make_fixture_db
