module sub_slice_number

! this file holds all the subroutines needed the main programs in this dir
! Qinya Liu, Caltech, May 2007

  implicit none

  include "constants.h"

contains


  subroutine tp2norm(ths,phs,thr,phr,thn,phn)

    ! find out the normal of the plane constructed
    ! by (ths,phs), (thr,phr), and the origin using
    ! cross-product.

    real(kind=CUSTOM_REAL) :: ths, phs, thr, phr, thn, phn

    real(kind=CUSTOM_REAL) :: xs, ys, zs, xr, yr, zr, nx, ny, nz

    call tp2xyz(ths,phs,xs,ys,zs)
    call tp2xyz(thr,phr,xr,yr,zr)
    nx = ys*zr - zs*yr
    ny = -(xs*zr - zs*xr)
    nz = xs*yr - ys*xr

    call xyz2tp(nx,ny,nz,thn,phn)


  end subroutine tp2norm

  ! ------------------------------------------------

  subroutine tp2norm2(ths,phs,thr,phr,nx,ny,nz)

! find out the normal of the plane (nx,ny,nz)
! constructed by (ths,phs), (thr,phr), and the origin.

    real(kind=CUSTOM_REAL) :: ths, phs, thr, phr, nx, ny, nz

    real(kind=CUSTOM_REAL) :: xs, ys, zs, xr, yr, zr, nr

    call tp2xyz(ths,phs,xs,ys,zs)
    call tp2xyz(thr,phr,xr,yr,zr)
    nx = ys*zr - zs*yr
    ny = -(xs*zr - zs*xr)
    nz = xs*yr - ys*xr
    nr = sqrt(nx*nx+ny*ny+nz*nz)
    nx = nx/nr
    ny = ny/nr
    nz = nz/nr

  end subroutine tp2norm2

  ! -------------------------------------------------

  subroutine norm_rot(thn,phn,th,ph,th_new,ph_new)


    ! coordinates change from (th,ph) to (th_new,ph_new)
    ! according to a rotation that converts (thn,phn) to
    ! z axis

    real(kind=CUSTOM_REAL) :: thn, phn, th, ph, th_new, ph_new

    real(kind=CUSTOM_REAL) :: rot(3,3), x, y, z, x_new, y_new, z_new

    rot(1,1) = cos(thn)*cos(phn)
    rot(1,2) = cos(thn)*sin(phn)
    rot(1,3) = -sin(thn)
    rot(2,1) = -sin(phn)
    rot(2,2) = cos(phn)
    rot(2,3) = 0
    rot(3,1) = sin(thn)*cos(phn)
    rot(3,2) = sin(thn)*sin(phn)
    rot(3,3) = cos(thn)

    call tp2xyz(th,ph,x,y,z)
    x_new = rot(1,1) * x + rot(1,2) * y + rot(1,3) * z
    y_new = rot(2,1) * x + rot(2,2) * y + rot(2,3) * z
    z_new = rot(3,1) * x + rot(3,2) * y + rot(3,3) * z

    call xyz2tp(x_new,y_new,z_new,th_new,ph_new)

  end subroutine norm_rot

  ! -------------------------------------------------

  subroutine linarray(array,n,a0,a1)


    real(kind=CUSTOM_REAL) :: array(*), a0, a1
    integer n

    real(kind=CUSTOM_REAL) :: da
    integer :: i

    if (n <= 0) stop 'number of points in linarray n > 0'
    if (n == 1 .and. abs(a0-a1) > EPS) stop 'number of points n = 1 asks for a0 = a1'

    if (n == 1) then
      array(1) = a0
    else
      da = (a1-a0)/(n-1)
      do i = 1, n
        array(i) = a0 + (i-1) * da
      enddo
    endif

  end subroutine linarray

  ! -------------------------------------------------

  subroutine rotmat(clon,clat,grot,rotation_matrix)

    ! this function calculate the 3x3 rotation matrix from the AB chunk
    ! frame to the actual frame defined by (clon,clat,grot)

    real(kind=CUSTOM_REAL) :: clon,clat,grot, rotation_matrix(3,3)

    real(kind=CUSTOM_REAL) :: alpha, beta, gamma, sina, cosa, sinb, cosb, &
               sing, cosg

    ! compute colatitude and longitude
    alpha = clon
    beta = pi/2 - clat
    gamma = grot

    sina = sin(alpha)
    cosa = cos(alpha)
    sinb = sin(beta)
    cosb = cos(beta)
    sing = sin(gamma)
    cosg = cos(gamma)

    ! define rotation matrix
    rotation_matrix(1,1) = cosg*cosb*cosa-sing*sina
    rotation_matrix(1,2) = -sing*cosb*cosa-cosg*sina
    rotation_matrix(1,3) = sinb*cosa
    rotation_matrix(2,1) = cosg*cosb*sina+sing*cosa
    rotation_matrix(2,2) = -sing*cosb*sina+cosg*cosa
    rotation_matrix(2,3) = sinb*sina
    rotation_matrix(3,1) = -cosg*sinb
    rotation_matrix(3,2) = sing*sinb
    rotation_matrix(3,3) = cosb

  end subroutine rotmat

  ! -------------------------------------------

  subroutine norm_rot_back(thn,phn,th,ph,th_old,ph_old)

    !coordinates change from (th_old,ph_old) to (th,ph)
    ! according to a rotation that converts (thn,phn) to
    ! z axis

    real(kind=CUSTOM_REAL) :: thn, phn, th, ph, th_old, ph_old

    real(kind=CUSTOM_REAL) :: rot(3,3), x,y,z, x_old, y_old, z_old

    rot(1,1) = cos(thn)*cos(phn)
    rot(1,2) = cos(thn)*sin(phn)
    rot(1,3) = -sin(thn)
    rot(2,1) = -sin(phn)
    rot(2,2) = cos(phn)
    rot(2,3) = 0
    rot(3,1) = sin(thn)*cos(phn)
    rot(3,2) = sin(thn)*sin(phn)
    rot(3,3) = cos(thn)
    call tp2xyz(th,ph,x,y,z)
    x_old = rot(1,1) * x + rot(2,1) * y + rot(3,1) * z
    y_old = rot(1,2) * x + rot(2,2) * y + rot(3,2) * z
    z_old = rot(1,3) * x + rot(2,3) * y + rot(3,3) * z

    call xyz2tp(x_old,y_old,z_old,th_old,ph_old)

  end subroutine norm_rot_back

  ! ----------------------------------------------

  subroutine tp2xyz(th,ph,x,y,z)

    real(CUSTOM_REAL) :: th,ph,x,y,z

    ! convert (th,ph) to (x,y,z) on unit sphere
    x = sin(th) * cos(ph)
    y = sin(th) * sin(ph)
    z = cos(th)

  end subroutine tp2xyz

  ! ---------------------------------------

  subroutine xyz2tp(x,y,z,th,ph)

    real(kind=CUSTOM_REAL) :: x, y, z, th, ph

    ! convert x,y,z to a point on unit sphere

    ph = atan2(y,x)
    th = atan2(sqrt(x*x+y*y),z)
    if (th < 0) th = th + pi
    if (ph < 0) ph = ph + 2*pi

  end subroutine xyz2tp

  ! ---------------------------------------

  subroutine chunk_map(k,xx,yy,zz,xi,eta)

    ! this program get the xi,eta for (xx,yy,zz)
    ! point under the k'th chunk coordinate
    ! transformation

    real(kind=CUSTOM_REAL) :: xx, yy, zz, xi, eta
    integer :: k

    real(kind=CUSTOM_REAL) :: x, y, z

    x = xx; y = yy; z = zz
    if (x < EPS .and. x >= 0) x = EPS
    if (x > -EPS .and. x < 0) x = -EPS
    if (y < EPS .and. y >= 0) y = EPS
    if (y > -EPS .and. y < 0) y = -EPS
    if (z < EPS .and. z >= 0) z = EPS
    if (z > -EPS .and. z < 0) z = -EPS

    if (k == 1) then ! CHUNK_AB
      xi = atan(y/z); eta = atan(-x/z)
      if (z < 0)  xi = 10
    else if (k == 2) then ! CHUNK_AC
      xi = atan(-z/y); eta = atan(x/y)
      if (y > 0) xi = 10
    else if (k == 3) then ! CHUNK_BC
      xi = atan(-z/x); eta = atan(-y/x)
      if (x > 0) xi = 10
    else if (k == 4) then ! CHUNK_AC'
      xi = atan(-z/y); eta = atan(-x/y)
      if (y < 0) xi = 10
    else if (k == 5) then ! CHUNK_BC'
      xi = atan(z/x); eta = atan(-y/x)
      if (x < 0) xi = 10
    else if (k == 6) then ! CHUNK_AB'
      xi = atan(y/z); eta = atan(x/z)
      if (z > 0)  xi = 10
    else
      stop 'chunk number k < 6'
    endif

  end subroutine chunk_map


  ! ------------------------------------------

  subroutine compact_int_array(s,np,ss,np2)

    integer :: s(*), ss(*)
    integer :: np, np2, i

    ss(1) = s(1)
    np2 = 1

    do i = 2, np
      if (all((s(i)-ss(1:np2)) /= 0)) then
        np2 = np2 + 1
        ss(np2) = s(i)
      endif
    enddo

  end subroutine compact_int_array


! -----------------------------------------


  subroutine distaz (sta,sto,epa,epo,dk,dd,daze,dazs)

    ! output: dk, dd -- distance in km and radians
    !         daze, dazs -- azimuth and back-azimuth

    real(kind=CUSTOM_REAL) :: sta, sto, epa, epo, dk, dd, daze, dazs

    real(kind=CUSTOM_REAL) :: rad, sa, ea, ssa, csa, so, eo, sea, cea, ces, &
               ses

    rad=pi/180.0d0

    sa  = atan(.993270*tan(sta*rad))
    ea  = atan(.993270*tan(epa*rad))
    ssa = sin(sa)
    csa = cos(sa)
    so  = sto*rad
    eo  = epo*rad
    sea = sin(ea)
    cea = cos(ea)
    ces = cos(eo-so)
    ses = sin(eo-so)

    if  (sa==ea) then
      if (sto==epo) then
        dk =0.00
        dd =0.00
        daze=0.0
        dazs=0.0
        return
      endif
    endif

    if (sta==90.) then
      if (epa==90.0) then
        dk =0.00
        dd =0.00
        daze=0.00
        dazs=0.00
        return
      endif
    endif

    if (sta==-90.0) then
      if (epa==-90.0) then
        dk =0.00
        dd =0.00
        daze=0.00
        dazs=0.00
        return
      endif
    endif

    dd = ssa*sea+csa*cea*ces
    if (dd /=0.) dd=atan(sqrt(1.0-dd*dd)/dd)
    if (dd==0.) dd=pi/2.0
    if (dd<0.0) dd=dd+pi
    dd = dd/rad
    dk = dd*111.19

    dazs = atan2(-ses,(ssa/csa*cea-sea*ces))
    daze = atan2(ses,(sea*csa/cea-ssa*ces))
    dazs = dazs/rad
    daze = daze/rad
    if (dazs<0.00)  dazs=dazs+360.0
    if (daze<0.00) daze=daze+360.0

  end subroutine distaz

! ---------------------------------------


  subroutine latlon2norm(lats,lons,latr,lonr,xs,ys,zs,nx,ny,nz,xc,yc,zc)

    ! [xs,ys,zs,nx,ny,nz] = latlon2norm(lats,lons,latr,lonr)
    ! find out the coordinates of the source (xs,ys,zs)
    ! the normal of the plane (nx,ny,nz)
    ! the center of the source and receiver garc (xc,yc,zc)
    ! constructed by (lons,lats), (lonr,latr), and the origin
    ! for a "banana cross-section"

    real(kind=CUSTOM_REAL) :: lats, lons, latr, lonr, xs, ys, zs, nx, ny, nz, xc, yc, zc

    real(kind=CUSTOM_REAL) :: ts, ps, tr, pr, xr, yr, zr, nr

    ts = (90 - lats) * pi/180
    ps = lons * pi / 180
    tr = (90 - latr) * pi/180
    pr = lonr * pi / 180

    call tp2xyz(ts,ps,xs,ys,zs)
    call tp2xyz(tr,pr,xr,yr,zr)

    nx = ys*zr - zs*yr
    ny = -(xs*zr - zs*xr)
    nz = xs*yr - ys*xr
    nr = sqrt(nx*nx+ny*ny+nz*nz)
    nx = nx/nr
    ny = ny/nr
    nz = nz/nr

    xc=(xr+xs)/2
    yc=(yr+ys)/2
    zc=(zr+zs)/2

  end subroutine latlon2norm

  ! ---------------------------------------




end module sub_slice_number
