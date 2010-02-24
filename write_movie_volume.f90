!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
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

! create file OUTPUT_FILES/values_from_mesher.h based upon DATA/Par_file
! in order to compile the solver with the right array sizes

!---------------------------------------------------------------------------------
! this subroutine counts the number of points and elements within the movie volume
! in this processor slice, and returns arrays that keep track of them, both in global and local indexing schemes

  subroutine count_points_movie_volume(prname,ibool_crust_mantle, xstore_crust_mantle,ystore_crust_mantle, &
                      zstore_crust_mantle,MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST,MOVIE_NORTH,MOVIE_SOUTH, &
                      MOVIE_COARSE,npoints_3dmovie,nspecel_3dmovie,num_ibool_3dmovie, &
                      mask_ibool_3dmovie,mask_3dmovie)

  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

! input
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle
  double precision :: MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST,MOVIE_NORTH,MOVIE_SOUTH
  logical :: MOVIE_COARSE
  character(len=150) :: prname

! output
  integer :: npoints_3dmovie,nspecel_3dmovie
  integer, dimension(NGLOB_CRUST_MANTLE) :: num_ibool_3dmovie
  logical, dimension(NGLOB_CRUST_MANTLE) :: mask_ibool_3dmovie
  logical, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: mask_3dmovie

! variables
  integer :: ipoints_3dmovie,ispecel_3dmovie,ispec,iglob,i,j,k,NIT
  real(kind=custom_real) :: rval,thetaval,phival

  if(MOVIE_COARSE) then
    NIT = NGLLX-1
  else
    NIT = 1
  endif
  ipoints_3dmovie=0
  num_ibool_3dmovie(:) = -99
  ispecel_3dmovie = 0
  mask_ibool_3dmovie(:)=.false.
  mask_3dmovie(:,:,:,:)=.false.
  ! create name of database
  open(unit=IOUT,file=trim(prname)//'movie3D_info.txt',status='unknown')

  !find and count points within given region for storing movie
      do ispec = 1,NSPEC_CRUST_MANTLE
        !output element if center of element is in the given region
        iglob    = ibool_crust_mantle((NGLLX+1)/2,(NGLLY+1)/2,(NGLLZ+1)/2,ispec)
        rval     = xstore_crust_mantle(iglob)
        thetaval = ystore_crust_mantle(iglob)
        phival   = zstore_crust_mantle(iglob)
      ! we alread changed xyz back to rthetaphi
        if( (rval < MOVIE_TOP .and. rval > MOVIE_BOTTOM) .and. &
            (thetaval > MOVIE_NORTH .and. thetaval < MOVIE_SOUTH) .and. &
            ( (phival < MOVIE_EAST .and. phival > MOVIE_WEST) .or. &
              ( (MOVIE_EAST < MOVIE_WEST) .and. (phival >MOVIE_EAST .or. phival < MOVIE_WEST) ) ) ) then
            ispecel_3dmovie=ispecel_3dmovie+1
              do k=1,NGLLZ,NIT
               do j=1,NGLLY,NIT
                do i=1,NGLLX,NIT
                 iglob    = ibool_crust_mantle(i,j,k,ispec)
                 if(.not. mask_ibool_3dmovie(iglob)) then
                  ipoints_3dmovie = ipoints_3dmovie + 1
                  mask_ibool_3dmovie(iglob)=.true.
                  mask_3dmovie(i,j,k,ispec)=.true.
                  num_ibool_3dmovie(iglob)= ipoints_3dmovie
                 endif
                enddo !i
               enddo !j
              enddo !k
        endif !in region
      enddo !ispec
   npoints_3dmovie=ipoints_3dmovie
   nspecel_3dmovie=ispecel_3dmovie

   write(IOUT,*) npoints_3dmovie, nspecel_3dmovie
   close(IOUT)

  end subroutine count_points_movie_volume

! -----------------------------------------------------------------
! writes meshfiles to merge with solver snapshots for 3D volume movies.  Also computes and outputs
! the rotation matrix nu_3dmovie required to transfer to a geographic coordinate system

  subroutine write_movie_volume_mesh(npoints_3dmovie,prname,ibool_crust_mantle,xstore_crust_mantle, &
                         ystore_crust_mantle,zstore_crust_mantle, muvstore_crust_mantle_3dmovie, &
                         mask_3dmovie,mask_ibool_3dmovie,num_ibool_3dmovie,nu_3dmovie,MOVIE_COARSE)

  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  !input
  integer :: npoints_3dmovie
  integer, dimension(NGLOB_CRUST_MANTLE) :: num_ibool_3dmovie
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: muvstore_crust_mantle_3dmovie
  character(len=150) :: prname
  logical, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: mask_3dmovie
  logical, dimension(NGLOB_CRUST_MANTLE) :: mask_ibool_3dmovie
  logical :: MOVIE_COARSE
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle

  !output
  real(kind=CUSTOM_REAL), dimension(3,3,npoints_3dmovie) :: nu_3dmovie
  real(kind=CUSTOM_REAL), dimension(npoints_3dmovie) :: store_val3D_mu

  !variables
  integer :: ipoints_3dmovie,ispecele,ispec,i,j,k,iglob,iglob1,iglob2,iglob3,iglob4,iglob5,iglob6,iglob7,iglob8
  integer :: n1,n2,n3,n4,n5,n6,n7,n8,NIT
  real(kind=CUSTOM_REAL) :: rval,thetaval,phival,xval,yval,zval,st,ct,sp,cp
  real(kind=CUSTOM_REAL), dimension(npoints_3dmovie) :: store_val3D_x,store_val3D_y, store_val3D_z

  if(NDIM /= 3) stop 'movie volume output requires NDIM = 3'

  if(MOVIE_COARSE) then
    NIT = NGLLX-1
  else
    NIT = 1
  endif

   ipoints_3dmovie=0
    do ispec=1,NSPEC_CRUST_MANTLE
     do k=1,NGLLZ,NIT
      do j=1,NGLLY,NIT
       do i=1,NGLLX,NIT
         if(mask_3dmovie(i,j,k,ispec)) then
          ipoints_3dmovie=ipoints_3dmovie+1
          iglob= ibool_crust_mantle(i,j,k,ispec)
          rval     = xstore_crust_mantle(iglob)
          thetaval = ystore_crust_mantle(iglob)
          phival   = zstore_crust_mantle(iglob)
!x,y,z store have been converted to r theta phi already, need to revert back for xyz output
          call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)
          store_val3D_x(ipoints_3dmovie)=xval
          store_val3D_y(ipoints_3dmovie)=yval
          store_val3D_z(ipoints_3dmovie)=zval
          store_val3D_mu(ipoints_3dmovie)=muvstore_crust_mantle_3dmovie(i,j,k,ispec)
          st = sin(thetaval)
          ct = cos(thetaval)
          sp = sin(phival)
          cp = cos(phival)
          nu_3dmovie(1,1,ipoints_3dmovie)=-ct*cp
          nu_3dmovie(1,2,ipoints_3dmovie)=-ct*sp
          nu_3dmovie(1,3,ipoints_3dmovie)=st
          nu_3dmovie(2,1,ipoints_3dmovie)=-sp
          nu_3dmovie(2,2,ipoints_3dmovie)=cp
          nu_3dmovie(2,3,ipoints_3dmovie)=0.d0
          nu_3dmovie(3,1,ipoints_3dmovie)=st*cp
          nu_3dmovie(3,2,ipoints_3dmovie)=st*sp
          nu_3dmovie(3,3,ipoints_3dmovie)=ct
        endif !mask_3dmovie
       enddo  !i
      enddo  !j
     enddo  !k
    enddo !ispec
   open(unit=IOUT,file=trim(prname)//'movie3D_x.bin',status='unknown',form='unformatted')
   if(npoints_3dmovie>0) then
     write(IOUT) store_val3D_x(1:npoints_3dmovie)
   endif
   close(IOUT)
   open(unit=IOUT,file=trim(prname)//'movie3D_y.bin',status='unknown',form='unformatted')
   if(npoints_3dmovie>0) then
     write(IOUT) store_val3D_y(1:npoints_3dmovie)
   endif
   close(IOUT)

   open(unit=IOUT,file=trim(prname)//'movie3D_z.bin',status='unknown',form='unformatted')
   if(npoints_3dmovie>0) then
     write(IOUT) store_val3D_z(1:npoints_3dmovie)
   endif
   close(IOUT)

   open(unit=IOUT,file=trim(prname)//'ascii_output.txt',status='unknown')
   if(npoints_3dmovie>0) then
     do i=1,npoints_3dmovie
       write(IOUT,*) store_val3D_x(i),store_val3D_y(i),store_val3D_z(i),store_val3D_mu(i)
     enddo
   endif
   close(IOUT)
   open(unit=IOUT,file=trim(prname)//'movie3D_elements.bin',status='unknown',form='unformatted')
   ispecele=0
 !  open(unit=IOUT,file=trim(prname)//'movie3D_elements.txt',status='unknown')
   do ispec=1,NSPEC_CRUST_MANTLE
    if(MOVIE_COARSE) then
      iglob=ibool_crust_mantle(1,1,1,ispec)
    else
      iglob=ibool_crust_mantle(3,3,3,ispec)
    endif
    if(mask_ibool_3dmovie(iglob)) then  !this element is in the region
     ispecele  = ispecele+1
     do k=1,NGLLZ-1,NIT
      do j=1,NGLLY-1,NIT
       do i=1,NGLLX-1,NIT
        ! if(mask_3dmovie(i,j,k,ispec)) then
          iglob1 = ibool_crust_mantle(i,j,k,ispec)
          iglob2 = ibool_crust_mantle(i+NIT,j,k,ispec)
          iglob3 = ibool_crust_mantle(i+NIT,j+NIT,k,ispec)
          iglob4 = ibool_crust_mantle(i,j+NIT,k,ispec)
          iglob5 = ibool_crust_mantle(i,j,k+NIT,ispec)
          iglob6 = ibool_crust_mantle(i+NIT,j,k+NIT,ispec)
          iglob7 = ibool_crust_mantle(i+NIT,j+NIT,k+NIT,ispec)
          iglob8 = ibool_crust_mantle(i,j+NIT,k+NIT,ispec)
          n1 = num_ibool_3dmovie(iglob1)-1
          n2 = num_ibool_3dmovie(iglob2)-1
          n3 = num_ibool_3dmovie(iglob3)-1
          n4 = num_ibool_3dmovie(iglob4)-1
          n5 = num_ibool_3dmovie(iglob5)-1
          n6 = num_ibool_3dmovie(iglob6)-1
          n7 = num_ibool_3dmovie(iglob7)-1
          n8 = num_ibool_3dmovie(iglob8)-1
          write(IOUT) n1,n2,n3,n4,n5,n6,n7,n8
        !  write(57,*) n1,n2,n3,n4,n5,n6,n7,n8
       !  endif !mask3dmovie
       enddo !i
      enddo !j
     enddo !k
    endif
    enddo !ispec
  close(IOUT)
!  close(57)
 end subroutine write_movie_volume_mesh

! ---------------------------------------------

  subroutine write_movie_volume_strains(myrank,npoints_3dmovie,LOCAL_PATH,MOVIE_VOLUME_TYPE,MOVIE_COARSE, &
                    it,eps_trace_over_3_crust_mantle,epsilondev_crust_mantle,muvstore_crust_mantle_3dmovie, &
                    mask_3dmovie,nu_3dmovie)


  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  ! input
  integer :: myrank,npoints_3dmovie,MOVIE_VOLUME_TYPE,it
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: eps_trace_over_3_crust_mantle
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: epsilondev_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: muvstore_crust_mantle_3dmovie
  logical, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: mask_3dmovie
  logical :: MOVIE_COARSE
  real(kind=CUSTOM_REAL), dimension(3,3,npoints_3dmovie) :: nu_3dmovie
  character(len=150) LOCAL_PATH,outputname

  ! variables
  character(len=150) prname
  integer :: ipoints_3dmovie,i,j,k,ispec,NIT
  real(kind=CUSTOM_REAL) :: muv_3dmovie
  real(kind=CUSTOM_REAL),dimension(3,3) :: eps_loc,eps_loc_new
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: store_val3d_NN,store_val3d_EE,store_val3d_ZZ,&
                                                     store_val3d_NE,store_val3d_NZ,store_val3d_EZ

  character(len=1) movie_prefix

  allocate(store_val3d_NN(npoints_3dmovie))
  allocate(store_val3d_EE(npoints_3dmovie))
  allocate(store_val3d_ZZ(npoints_3dmovie))
  allocate(store_val3d_NE(npoints_3dmovie))
  allocate(store_val3d_NZ(npoints_3dmovie))
  allocate(store_val3d_EZ(npoints_3dmovie))

  if(NDIM /= 3) call exit_MPI(myrank, 'write_movie_volume requires NDIM = 3')

  if(MOVIE_VOLUME_TYPE == 1) then
      movie_prefix='E' ! strain
  else if(MOVIE_VOLUME_TYPE == 2) then
      movie_prefix='S' ! time integral of strain
  else if(MOVIE_VOLUME_TYPE == 3) then
      movie_prefix='P' ! potency, or integral of strain x \mu
  endif
  if(MOVIE_COARSE) then
   NIT = NGLLX-1
  else
   NIT = 1
  endif
  
  write(prname,"('proc',i6.6)") myrank
  ipoints_3dmovie=0
  do ispec=1,NSPEC_CRUST_MANTLE
   do k=1,NGLLZ,NIT
    do j=1,NGLLY,NIT
     do i=1,NGLLX,NIT
      if(mask_3dmovie(i,j,k,ispec)) then
       ipoints_3dmovie=ipoints_3dmovie+1
       muv_3dmovie=muvstore_crust_mantle_3dmovie(i,j,k,ispec)
       eps_loc(1,1)=eps_trace_over_3_crust_mantle(i,j,k,ispec) + epsilondev_crust_mantle(1,i,j,k,ispec)
       eps_loc(2,2)=eps_trace_over_3_crust_mantle(i,j,k,ispec) + epsilondev_crust_mantle(2,i,j,k,ispec)
       eps_loc(3,3)=eps_trace_over_3_crust_mantle(i,j,k,ispec)- &
                 epsilondev_crust_mantle(1,i,j,k,ispec) - epsilondev_crust_mantle(2,i,j,k,ispec)
       eps_loc(1,2)=epsilondev_crust_mantle(3,i,j,k,ispec)
       eps_loc(1,3)=epsilondev_crust_mantle(4,i,j,k,ispec)
       eps_loc(2,3)=epsilondev_crust_mantle(5,i,j,k,ispec)
       eps_loc(2,1)=eps_loc(1,2)
       eps_loc(3,1)=eps_loc(1,3)
       eps_loc(3,2)=eps_loc(2,3)

  ! rotate eps_loc to spherical coordinates
    eps_loc_new(:,:) = matmul(matmul(nu_3dmovie(:,:,ipoints_3dmovie),eps_loc(:,:)), transpose(nu_3dmovie(:,:,ipoints_3dmovie)))
       if(MOVIE_VOLUME_TYPE == 3) eps_loc_new(:,:) = eps_loc(:,:)*muv_3dmovie
       store_val3d_NN(ipoints_3dmovie)=eps_loc_new(1,1)
       store_val3d_EE(ipoints_3dmovie)=eps_loc_new(2,2)
       store_val3d_ZZ(ipoints_3dmovie)=eps_loc_new(3,3)
       store_val3d_NE(ipoints_3dmovie)=eps_loc_new(1,2)
       store_val3d_NZ(ipoints_3dmovie)=eps_loc_new(1,3)
       store_val3d_EZ(ipoints_3dmovie)=eps_loc_new(2,3)
      endif
     enddo
    enddo
   enddo
  enddo
  if(ipoints_3dmovie /= npoints_3dmovie) stop 'did not find the right number of points for 3D movie'

  write(outputname,"('proc',i6.6,'_movie3D_',a,'NN',i6.6,'.bin')") myrank,movie_prefix,it
  open(unit=27,file=trim(LOCAL_PATH)//'/'//trim(outputname),status='unknown',form='unformatted')
  write(27) store_val3d_NN(1:npoints_3dmovie)
  close(27)

  write(outputname,"('proc',i6.6,'_movie3D_',a,'EE',i6.6,'.bin')") myrank,movie_prefix,it
  open(unit=27,file=trim(LOCAL_PATH)//'/'//trim(outputname),status='unknown',form='unformatted')
  write(27) store_val3d_EE(1:npoints_3dmovie)
  close(27)

  write(outputname,"('proc',i6.6,'_movie3D_',a,'ZZ',i6.6,'.bin')") myrank,movie_prefix,it
  open(unit=27,file=trim(LOCAL_PATH)//'/'//trim(outputname),status='unknown',form='unformatted')
  write(27) store_val3d_ZZ(1:npoints_3dmovie)
  close(27)

  write(outputname,"('proc',i6.6,'_movie3D_',a,'NE',i6.6,'.bin')") myrank,movie_prefix,it
  open(unit=27,file=trim(LOCAL_PATH)//'/'//trim(outputname),status='unknown',form='unformatted')
  write(27) store_val3d_NE(1:npoints_3dmovie)
  close(27)


  write(outputname,"('proc',i6.6,'_movie3D_',a,'NZ',i6.6,'.bin')") myrank,movie_prefix,it
  open(unit=27,file=trim(LOCAL_PATH)//'/'//trim(outputname),status='unknown',form='unformatted')
  write(27) store_val3d_NZ(1:npoints_3dmovie)
  close(27)

  write(outputname,"('proc',i6.6,'_movie3D_',a,'EZ',i6.6,'.bin')") myrank,movie_prefix,it
  open(unit=27,file=trim(LOCAL_PATH)//'/'//trim(outputname),status='unknown',form='unformatted')
  write(27) store_val3d_EZ(1:npoints_3dmovie)
  close(27)

  end subroutine write_movie_volume_strains

! ---------------------------------------------
  subroutine write_movie_volume_vector(myrank,it,npoints_3dmovie,LOCAL_PATH,MOVIE_VOLUME_TYPE, &
                    MOVIE_COARSE,ibool_crust_mantle,vector_crust_mantle,scalingval,mask_3dmovie,nu_3dmovie)
  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  ! input
  integer :: myrank,npoints_3dmovie,MOVIE_VOLUME_TYPE,it
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  real(kind=CUSTOM_REAL), dimension(3,NGLOB_CRUST_MANTLE) :: vector_crust_mantle,vector_scaled
  real(kind=CUSTOM_REAL), dimension(3,3,npoints_3dmovie) :: nu_3dmovie
  double precision :: scalingval
  real(kind=CUSTOM_REAL), dimension(3) :: vector_local,vector_local_new
  logical, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: mask_3dmovie
  logical :: MOVIE_COARSE
  character(len=150) LOCAL_PATH

  ! variables
  integer :: ipoints_3dmovie,i,j,k,ispec,NIT,iglob
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: store_val3d_N,store_val3d_E,store_val3d_Z
  character(len=150) outputname
  character(len=2) movie_prefix

  if(NDIM /= 3) call exit_MPI(myrank,'write_movie_volume requires NDIM = 3')

  allocate(store_val3d_N(npoints_3dmovie))
  allocate(store_val3d_E(npoints_3dmovie))
  allocate(store_val3d_Z(npoints_3dmovie))

  if(MOVIE_VOLUME_TYPE == 5) then
      movie_prefix='DI' ! displacement
  else if(MOVIE_VOLUME_TYPE == 6) then
      movie_prefix='VE' ! velocity
  endif
  if(MOVIE_COARSE) then
   NIT = NGLLX-1
  else
   NIT = 1
  endif

  if(CUSTOM_REAL == SIZE_REAL) then
    vector_scaled = vector_crust_mantle*sngl(scalingval)
  else
    vector_scaled = vector_crust_mantle*scalingval
  endif

  ipoints_3dmovie=0
  do ispec=1,NSPEC_CRUST_MANTLE
   do k=1,NGLLZ,NIT
    do j=1,NGLLY,NIT
     do i=1,NGLLX,NIT
      if(mask_3dmovie(i,j,k,ispec)) then
       ipoints_3dmovie=ipoints_3dmovie+1
       iglob = ibool_crust_mantle(i,j,k,ispec)
       vector_local(:) = vector_scaled(:,iglob)

  ! rotate eps_loc to spherical coordinates
       vector_local_new(:) = matmul(nu_3dmovie(:,:,ipoints_3dmovie), vector_local(:))
       store_val3d_N(ipoints_3dmovie)=vector_local_new(1)
       store_val3d_E(ipoints_3dmovie)=vector_local_new(2)
       store_val3d_Z(ipoints_3dmovie)=vector_local_new(3)
      endif
     enddo
    enddo
   enddo
  enddo
  close(IOUT)
  if(ipoints_3dmovie /= npoints_3dmovie) stop 'did not find the right number of points for 3D movie'

  write(outputname,"('proc',i6.6,'_movie3D_',a,'N',i6.6,'.bin')") myrank,movie_prefix,it
  open(unit=27,file=trim(LOCAL_PATH)//'/'//trim(outputname),status='unknown',form='unformatted')
  write(27) store_val3d_N(1:npoints_3dmovie)
  close(27)

  write(outputname,"('proc',i6.6,'_movie3D_',a,'E',i6.6,'.bin')") myrank,movie_prefix,it
  open(unit=27,file=trim(LOCAL_PATH)//'/'//trim(outputname),status='unknown',form='unformatted')
  write(27) store_val3d_E(1:npoints_3dmovie)
  close(27)

  write(outputname,"('proc',i6.6,'_movie3D_',a,'Z',i6.6,'.bin')") myrank,movie_prefix,it
  open(unit=27,file=trim(LOCAL_PATH)//'/'//trim(outputname),status='unknown',form='unformatted')
  write(27) store_val3d_Z(1:npoints_3dmovie)
  close(27)


  end subroutine write_movie_volume_vector

!--------------------

 subroutine write_movie_volume_divcurl(myrank,it,eps_trace_over_3_crust_mantle,&
          accel_outer_core,kappavstore_outer_core,rhostore_outer_core,ibool_outer_core, &
          eps_trace_over_3_inner_core,epsilondev_crust_mantle,&
          epsilondev_inner_core)
    include "constants.h"
    include "OUTPUT_FILES/values_from_mesher.h"
    ! div
    integer :: myrank,it,ispec,iglob,i,j,k
    real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_STR_OR_ATT) :: epsilondev_inner_core
    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_STRAIN_ONLY) :: eps_trace_over_3_inner_core
    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ADJOINT) :: rhostore_outer_core, &
                            kappavstore_outer_core,ibool_outer_core
    real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: accel_outer_core
    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: eps_trace_over_3_crust_mantle
    real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: epsilondev_crust_mantle
    real(kind=CUSTOM_REAL) :: rhol,kappal
    real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: div_s_outer_core

    character(len=150) LOCAL_PATH,outputname


    write(outputname,"('proc',i6.6,'_crust_mantle_div_displ_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//trim(outputname),status='unknown',form='unformatted')
    write(27) eps_trace_over_3_crust_mantle
    close(27)

! we use div s = - p / kappa = rhostore_outer_core * accel_outer_core / kappavstore_outer_core
    allocate(div_s_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ADJOINT))
    do ispec = 1, NSPEC_OUTER_CORE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_outer_core(i,j,k,ispec)
            rhol = rhostore_outer_core(i,j,k,ispec)
            kappal = kappavstore_outer_core(i,j,k,ispec)
            div_s_outer_core(i,j,k,ispec) = rhol * accel_outer_core(iglob) / kappal
          enddo
        enddo
      enddo
    enddo

    write(outputname,"('proc',i6.6,'_outer_core_div_displ_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//trim(outputname),status='unknown',form='unformatted')
    write(27)  div_s_outer_core
    close(27)

    deallocate(div_s_outer_core)


  !  write(outputname,"('proc',i6.6,'_outer_core_div_displ_it',i6.6,'.bin')") myrank,it
  !  open(unit=27,file=trim(LOCAL_PATH)//trim(outputname),status='unknown',form='unformatted')
  !  write(27)  ONE_THIRD * div_displ_outer_core
  !  close(27)

    write(outputname,"('proc',i6.6,'_inner_core_div_displ_proc_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//trim(outputname),status='unknown',form='unformatted')
    write(27) eps_trace_over_3_inner_core
    close(27)

! epsilondev

    write(outputname,"('proc',i6.6,'_crust_mantle_epsdev_displ_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//trim(outputname),status='unknown',form='unformatted')
    write(27) epsilondev_crust_mantle
    close(27)

    write(outputname,"('proc',i6.6,'inner_core_epsdev_displ_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//trim(outputname),status='unknown',form='unformatted')
    write(27) epsilondev_inner_core
    close(27)


  end subroutine write_movie_volume_divcurl


