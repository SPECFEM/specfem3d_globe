!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 5
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology July 2004
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine get_flags_boundaries(myrank,iregion_code,nspec,iproc_xi,iproc_eta,ispec, &
             xstore,ystore,zstore,iboun,iMPIcut_xi,iMPIcut_eta,ichunk, &
             idoubling,NPROC_XI,NPROC_ETA,rotation_matrix,ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD)

  implicit none

  include "constants.h"

  integer nspec
  integer ispec,ichunk,iregion_code,myrank,idoubling
  integer NPROC_XI,NPROC_ETA

  logical iboun(6,nspec)
  logical iMPIcut_xi(2,nspec),iMPIcut_eta(2,nspec)

  double precision xstore(NGLLX,NGLLY,NGLLZ)
  double precision ystore(NGLLX,NGLLY,NGLLZ)
  double precision zstore(NGLLX,NGLLY,NGLLZ)

  double precision ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD

  integer ia

! use iproc_xi and iproc_eta to determine MPI cut planes along xi and eta
  integer iproc_xi,iproc_eta

  double precision target,sizeslice
  double precision xelm(8),yelm(8),zelm(8)
  double precision xi(8),eta(8)
  double precision r1,r2,r3,r4,r5,r6,r7,r8
  double precision val_divide,radius_cube

! rotation matrix from Euler angles
  integer i,j
  double precision rotation_matrix(3,3)
  double precision vector_ori(3),vector_rotated(3)

! find the coordinates of the eight corner nodes of the element
  xelm(1)=xstore(1,1,1)
  yelm(1)=ystore(1,1,1)
  zelm(1)=zstore(1,1,1)
  xelm(2)=xstore(NGLLX,1,1)
  yelm(2)=ystore(NGLLX,1,1)
  zelm(2)=zstore(NGLLX,1,1)
  xelm(3)=xstore(NGLLX,NGLLY,1)
  yelm(3)=ystore(NGLLX,NGLLY,1)
  zelm(3)=zstore(NGLLX,NGLLY,1)
  xelm(4)=xstore(1,NGLLY,1)
  yelm(4)=ystore(1,NGLLY,1)
  zelm(4)=zstore(1,NGLLY,1)
  xelm(5)=xstore(1,1,NGLLZ)
  yelm(5)=ystore(1,1,NGLLZ)
  zelm(5)=zstore(1,1,NGLLZ)
  xelm(6)=xstore(NGLLX,1,NGLLZ)
  yelm(6)=ystore(NGLLX,1,NGLLZ)
  zelm(6)=zstore(NGLLX,1,NGLLZ)
  xelm(7)=xstore(NGLLX,NGLLY,NGLLZ)
  yelm(7)=ystore(NGLLX,NGLLY,NGLLZ)
  zelm(7)=zstore(NGLLX,NGLLY,NGLLZ)
  xelm(8)=xstore(1,NGLLY,NGLLZ)
  yelm(8)=ystore(1,NGLLY,NGLLZ)
  zelm(8)=zstore(1,NGLLY,NGLLZ)

! rotate origin of chunk back to North pole using transpose of rotation matrix
  if(NCHUNKS /= 6) then
  do ia=1,8
    vector_ori(1) = xelm(ia)
    vector_ori(2) = yelm(ia)
    vector_ori(3) = zelm(ia)
    do i=1,3
      vector_rotated(i) = ZERO
      do j=1,3
        vector_rotated(i)=vector_rotated(i)+rotation_matrix(j,i)*vector_ori(j)
      enddo
    enddo
    xelm(ia) = vector_rotated(1)
    yelm(ia) = vector_rotated(2)
    zelm(ia) = vector_rotated(3)
  enddo
  endif

! convert the 6 chunks to the cubed sphere coordinates xi and eta

! loop on the 8 control points
  do ia=1,8

! map back to cubed sphere, making sure we never divide by zero
  if((idoubling /= IFLAG_IN_CENTRAL_CUBE .and. &
     idoubling /= IFLAG_BOTTOM_CENTRAL_CUBE .and. &
     idoubling /= IFLAG_TOP_CENTRAL_CUBE .and. &
     idoubling /= IFLAG_IN_FICTITIOUS_CUBE) .or. (iregion_code == IREGION_INNER_CORE .and. INFLATE_CENTRAL_CUBE)) then

    if(ichunk == CHUNK_AB) then
      val_divide = zelm(ia)
      if(val_divide>-SMALL_VAL_ANGLE .and. val_divide <= ZERO) val_divide=-SMALL_VAL_ANGLE
      if(val_divide<SMALL_VAL_ANGLE .and. val_divide >= ZERO) val_divide=SMALL_VAL_ANGLE
      xi(ia)=datan(yelm(ia)/val_divide)
      eta(ia)=datan(-xelm(ia)/val_divide)

    else if(ichunk == CHUNK_AB_ANTIPODE) then
      val_divide = zelm(ia)
      if(val_divide>-SMALL_VAL_ANGLE .and. val_divide <= ZERO) val_divide=-SMALL_VAL_ANGLE
      if(val_divide<SMALL_VAL_ANGLE .and. val_divide >= ZERO) val_divide=SMALL_VAL_ANGLE
      xi(ia)=datan(yelm(ia)/val_divide)
      eta(ia)=datan(xelm(ia)/val_divide)

    else if(ichunk == CHUNK_AC) then
      val_divide = yelm(ia)
      if(val_divide>-SMALL_VAL_ANGLE .and. val_divide <= ZERO) val_divide=-SMALL_VAL_ANGLE
      if(val_divide<SMALL_VAL_ANGLE .and. val_divide >= ZERO) val_divide=SMALL_VAL_ANGLE
      xi(ia)=datan(-zelm(ia)/val_divide)
      eta(ia)=datan(xelm(ia)/val_divide)

    else if(ichunk == CHUNK_AC_ANTIPODE) then
      val_divide = yelm(ia)
      if(val_divide>-SMALL_VAL_ANGLE .and. val_divide <= ZERO) val_divide=-SMALL_VAL_ANGLE
      if(val_divide<SMALL_VAL_ANGLE .and. val_divide >= ZERO) val_divide=SMALL_VAL_ANGLE
      xi(ia)=datan(-zelm(ia)/val_divide)
      eta(ia)=datan(-xelm(ia)/val_divide)

    else if(ichunk == CHUNK_BC) then
      val_divide = xelm(ia)
      if(val_divide>-SMALL_VAL_ANGLE .and. val_divide <= ZERO) val_divide=-SMALL_VAL_ANGLE
      if(val_divide<SMALL_VAL_ANGLE .and. val_divide >= ZERO) val_divide=SMALL_VAL_ANGLE
      xi(ia)=datan(-zelm(ia)/val_divide)
      eta(ia)=datan(-yelm(ia)/val_divide)

    else if(ichunk == CHUNK_BC_ANTIPODE) then
      val_divide = xelm(ia)
      if(val_divide>-SMALL_VAL_ANGLE .and. val_divide <= ZERO) val_divide=-SMALL_VAL_ANGLE
      if(val_divide<SMALL_VAL_ANGLE .and. val_divide >= ZERO) val_divide=SMALL_VAL_ANGLE
      xi(ia)=datan(zelm(ia)/val_divide)
      eta(ia)=datan(-yelm(ia)/val_divide)
    else
      call exit_MPI(myrank,'incorrect chunk numbering in cubed sphere')
    endif

! map back by inverting central cube mapping if element inside the cube
  else

! size of the cube along cartesian axes before rotation
    radius_cube = R_CENTRAL_CUBE / dsqrt(3.d0)

    if(ichunk == CHUNK_AB) then
      xi(ia)=datan(yelm(ia)/radius_cube)
      eta(ia)=datan(-xelm(ia)/radius_cube)

    else if(ichunk == CHUNK_AB_ANTIPODE) then
      xi(ia)=datan(-yelm(ia)/radius_cube)
      eta(ia)=datan(-xelm(ia)/radius_cube)

    else if(ichunk == CHUNK_AC) then
      xi(ia)=datan(zelm(ia)/radius_cube)
      eta(ia)=datan(-xelm(ia)/radius_cube)

    else if(ichunk == CHUNK_AC_ANTIPODE) then
      xi(ia)=datan(-zelm(ia)/radius_cube)
      eta(ia)=datan(-xelm(ia)/radius_cube)

    else if(ichunk == CHUNK_BC) then
      xi(ia)=datan(zelm(ia)/radius_cube)
      eta(ia)=datan(yelm(ia)/radius_cube)

    else if(ichunk == CHUNK_BC_ANTIPODE) then
      xi(ia)=datan(zelm(ia)/radius_cube)
      eta(ia)=datan(-yelm(ia)/radius_cube)

    else
      call exit_MPI(myrank,'incorrect chunk numbering in central cube')
    endif

  endif

  enddo

  r1=dsqrt(xelm(1)**2+yelm(1)**2+zelm(1)**2)
  r2=dsqrt(xelm(2)**2+yelm(2)**2+zelm(2)**2)
  r3=dsqrt(xelm(3)**2+yelm(3)**2+zelm(3)**2)
  r4=dsqrt(xelm(4)**2+yelm(4)**2+zelm(4)**2)
  r5=dsqrt(xelm(5)**2+yelm(5)**2+zelm(5)**2)
  r6=dsqrt(xelm(6)**2+yelm(6)**2+zelm(6)**2)
  r7=dsqrt(xelm(7)**2+yelm(7)**2+zelm(7)**2)
  r8=dsqrt(xelm(8)**2+yelm(8)**2+zelm(8)**2)

! ****************************************************
!     determine if the element falls on a boundary
! ****************************************************

  iboun(:,ispec)=.false.

  target= (ANGULAR_WIDTH_XI_RAD/2.)*(ONE-SMALLVAL)

! on boundary 1: x=xmin
  if(xi(1)<-target .and. xi(4)<-target .and. xi(5)<-target .and. xi(8)<-target) &
    iboun(1,ispec)=.true.

! on boundary 2: xmax
  if(xi(2)>target .and. xi(3)>target .and. xi(6)>target .and. xi(7)>target) &
    iboun(2,ispec)=.true.

  target= (ANGULAR_WIDTH_ETA_RAD/2.)*(ONE-SMALLVAL)

! on boundary 3: ymin
  if(eta(1)<-target .and. eta(2)<-target .and. eta(5)<-target .and. eta(6)<-target) &
    iboun(3,ispec)=.true.

! on boundary 4: ymax
  if(eta(3)>target .and. eta(4)>target .and. eta(7)>target .and. eta(8)>target) &
    iboun(4,ispec)=.true.

! on boundary 5: bottom

  if(iregion_code /= IREGION_INNER_CORE .or. (iregion_code == IREGION_INNER_CORE .and. INFLATE_CENTRAL_CUBE)) then

    select case(iregion_code)
      case(IREGION_CRUST_MANTLE)
        target = (RCMB/R_EARTH)*(ONE + SMALLVAL)
      case(IREGION_OUTER_CORE)
        target = (RICB/R_EARTH)*(ONE + SMALLVAL)
      case(IREGION_INNER_CORE)
        target = R_CENTRAL_CUBE*(ONE + SMALLVAL)
    end select

    if(r1<target .and. r2<target .and. r3<target .and. r4<target) &
      iboun(5,ispec)=.true.

  else

! detect flat central cube boundary at the bottom instead of cubed sphere
! exclude elements that are located in the central cube itself

  if(idoubling /= IFLAG_IN_CENTRAL_CUBE .and. &
     idoubling /= IFLAG_BOTTOM_CENTRAL_CUBE .and. &
     idoubling /= IFLAG_TOP_CENTRAL_CUBE .and. &
     idoubling /= IFLAG_IN_FICTITIOUS_CUBE) then

  select case(ichunk)
    case(CHUNK_AB)
      target = (R_CENTRAL_CUBE/dsqrt(3.d0))*(ONE + SMALLVAL)
      if(zelm(1) < target .and. zelm(2) < target .and. zelm(3) < target .and. zelm(4) < target) iboun(5,ispec)=.true.

    case(CHUNK_AC)
      target = - (R_CENTRAL_CUBE/dsqrt(3.d0))*(ONE + SMALLVAL)
      if(yelm(1) > target .and. yelm(2) > target .and. yelm(3) > target .and. yelm(4) > target) iboun(5,ispec)=.true.

    case(CHUNK_BC)
      target = - (R_CENTRAL_CUBE/dsqrt(3.d0))*(ONE + SMALLVAL)
      if(xelm(1) > target .and. xelm(2) > target .and. xelm(3) > target .and. xelm(4) > target) iboun(5,ispec)=.true.

    case(CHUNK_AB_ANTIPODE)
      target = - (R_CENTRAL_CUBE/dsqrt(3.d0))*(ONE + SMALLVAL)
      if(zelm(1) > target .and. zelm(2) > target .and. zelm(3) > target .and. zelm(4) > target) iboun(5,ispec)=.true.

    case(CHUNK_AC_ANTIPODE)
      target = (R_CENTRAL_CUBE/dsqrt(3.d0))*(ONE + SMALLVAL)
      if(yelm(1) < target .and. yelm(2) < target .and. yelm(3) < target .and. yelm(4) < target) iboun(5,ispec)=.true.

    case(CHUNK_BC_ANTIPODE)
      target = (R_CENTRAL_CUBE/dsqrt(3.d0))*(ONE + SMALLVAL)
      if(xelm(1) < target .and. xelm(2) < target .and. xelm(3) < target .and. xelm(4) < target) iboun(5,ispec)=.true.

    case default
      call exit_MPI(myrank,'incorrect chunk code in flag detection')
  end select

  endif

  endif

! on boundary 6: top
  select case(iregion_code)
    case(IREGION_CRUST_MANTLE)
      target = R_UNIT_SPHERE*(ONE - SMALLVAL)
    case(IREGION_OUTER_CORE)
      target = (RCMB/R_EARTH)*(ONE - SMALLVAL)
    case(IREGION_INNER_CORE)
      target = (RICB/R_EARTH)*(ONE - SMALLVAL)
    case default
      call exit_MPI(myrank,'incorrect region code')
  end select
  if(r5>target .and. r6>target .and. r7>target .and. r8>target) &
    iboun(6,ispec)=.true.

! *******************************************************************
!     determine if the element falls on an MPI cut plane along xi
! *******************************************************************

! detect the MPI cut planes along xi in the cubed sphere

  iMPIcut_xi(:,ispec)=.false.

! angular size of a slice along xi
  sizeslice = ANGULAR_WIDTH_XI_RAD / NPROC_XI

! left cut-plane in the current slice along X = constant (Xmin of this slice)

  target = - (ANGULAR_WIDTH_XI_RAD/2.) + iproc_xi*sizeslice

! add geometrical tolerance
  target = target + SMALLVAL

  if(xi(1)<target .and. xi(4)<target .and. xi(5)<target .and. xi(8)<target) &
    iMPIcut_xi(1,ispec)=.true.

!  right cut-plane in the current slice along X = constant (Xmax of this slice)

  target = - (ANGULAR_WIDTH_XI_RAD/2.) + (iproc_xi+1)*sizeslice

! add geometrical tolerance
  target = target - SMALLVAL

  if(xi(2)>target .and. xi(3)>target .and. xi(6)>target .and. xi(7)>target) &
    iMPIcut_xi(2,ispec)=.true.

! ********************************************************************
!     determine if the element falls on an MPI cut plane along eta
! ********************************************************************

  iMPIcut_eta(:,ispec)=.false.

! angular size of a slice along eta
  sizeslice = ANGULAR_WIDTH_ETA_RAD / NPROC_ETA

! left cut-plane in the current slice along Y = constant (Ymin of this slice)

  target = - (ANGULAR_WIDTH_ETA_RAD/2.) + iproc_eta*sizeslice

! add geometrical tolerance
  target = target + SMALLVAL

  if(eta(1)<target .and. eta(2)<target .and. eta(5)<target .and. eta(6)<target) &
    iMPIcut_eta(1,ispec)=.true.

! right cut-plane in the current slice along Y = constant (Ymax of this slice)

  target = - (ANGULAR_WIDTH_ETA_RAD/2.) + (iproc_eta+1)*sizeslice

! add geometrical tolerance
  target = target - SMALLVAL

  if(eta(3)>target .and. eta(4)>target .and. eta(7)>target .and. eta(8)>target) &
    iMPIcut_eta(2,ispec)=.true.

  end subroutine get_flags_boundaries

