!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
!                            April 2011
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

! this subroutine counts the number of points and elements within the movie volume
! in this processor slice, and returns arrays that keep track of them, both in global and local indexing schemes

  subroutine count_points_movie_volume(prname,ibool_crust_mantle, xstore_crust_mantle,ystore_crust_mantle, &
                      zstore_crust_mantle,MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST,MOVIE_NORTH,MOVIE_SOUTH, &
                      MOVIE_COARSE,npoints_3dmovie,nspecel_3dmovie,num_ibool_3dmovie, &
                      mask_ibool,mask_3dmovie,myrank)

  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

! input
  integer :: myrank
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle
  double precision :: MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST,MOVIE_NORTH,MOVIE_SOUTH
  logical :: MOVIE_COARSE
  character(len=150) :: prname

! output
  integer :: npoints_3dmovie,nspecel_3dmovie
  integer, dimension(NGLOB_CRUST_MANTLE) :: num_ibool_3dmovie
  logical, dimension(NGLOB_CRUST_MANTLE) :: mask_ibool
  logical, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: mask_3dmovie

  ! local parameters
  integer :: ipoints_3dmovie,ispecel_3dmovie,ispec,iglob,i,j,k,iNIT,ier
  real(kind=custom_real) :: rval,thetaval,phival

  if(MOVIE_COARSE) then
    iNIT = NGLLX-1
  else
    iNIT = 1
  endif

  ipoints_3dmovie=0
  num_ibool_3dmovie(:) = -99
  ispecel_3dmovie = 0
  mask_ibool(:)=.false.
  mask_3dmovie(:,:,:,:)=.false.

  ! create name of database
  open(unit=IOUT,file=trim(prname)//'movie3D_info.txt', &
        status='unknown',iostat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error opening file movie3D_info.txt')

  !find and count points within given region for storing movie
  do ispec = 1,NSPEC_CRUST_MANTLE
    !output element if center of element is in the given region
    iglob    = ibool_crust_mantle((NGLLX+1)/2,(NGLLY+1)/2,(NGLLZ+1)/2,ispec)
    rval     = xstore_crust_mantle(iglob)
    thetaval = ystore_crust_mantle(iglob)
    phival   = zstore_crust_mantle(iglob)

    ! we already changed xyz back to rthetaphi
    if( (rval < MOVIE_TOP .and. rval > MOVIE_BOTTOM) .and. &
       (thetaval > MOVIE_NORTH .and. thetaval < MOVIE_SOUTH) .and. &
       ( (phival < MOVIE_EAST .and. phival > MOVIE_WEST) .or. &
       ( (MOVIE_EAST < MOVIE_WEST) .and. (phival >MOVIE_EAST .or. phival < MOVIE_WEST) ) ) ) then
      ispecel_3dmovie=ispecel_3dmovie+1
      do k=1,NGLLZ,iNIT
        do j=1,NGLLY,iNIT
          do i=1,NGLLX,iNIT
            iglob    = ibool_crust_mantle(i,j,k,ispec)
            if(.not. mask_ibool(iglob)) then
              ipoints_3dmovie = ipoints_3dmovie + 1
              mask_ibool(iglob)=.true.
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

!
!-------------------------------------------------------------------------------------------------
!

! writes meshfiles to merge with solver snapshots for 3D volume movies.  Also computes and outputs
! the rotation matrix nu_3dmovie required to transfer to a geographic coordinate system

  subroutine write_movie_volume_mesh(npoints_3dmovie,prname,ibool_crust_mantle,xstore_crust_mantle, &
                         ystore_crust_mantle,zstore_crust_mantle, muvstore_crust_mantle_3dmovie, &
                         mask_3dmovie,mask_ibool,num_ibool_3dmovie,nu_3dmovie,MOVIE_COARSE)

  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  !input
  integer :: npoints_3dmovie
  integer, dimension(NGLOB_CRUST_MANTLE) :: num_ibool_3dmovie
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_3DMOVIE) :: muvstore_crust_mantle_3dmovie
  character(len=150) :: prname
  logical, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: mask_3dmovie
  logical, dimension(NGLOB_CRUST_MANTLE) :: mask_ibool
  logical :: MOVIE_COARSE
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle

  !output
  real(kind=CUSTOM_REAL), dimension(3,3,npoints_3dmovie) :: nu_3dmovie
  real(kind=CUSTOM_REAL), dimension(npoints_3dmovie) :: store_val3D_mu

  !variables
  integer :: ipoints_3dmovie,ispecele,ispec,i,j,k,iglob,iglob1,iglob2,iglob3,iglob4,iglob5,iglob6,iglob7,iglob8
  integer :: n1,n2,n3,n4,n5,n6,n7,n8,iNIT
  real(kind=CUSTOM_REAL) :: rval,thetaval,phival,xval,yval,zval,st,ct,sp,cp
  real(kind=CUSTOM_REAL), dimension(npoints_3dmovie) :: store_val3D_x,store_val3D_y, store_val3D_z

  if(NDIM /= 3) stop 'movie volume output requires NDIM = 3'

  if(MOVIE_COARSE) then
    iNIT = NGLLX-1
  else
    iNIT = 1
  endif

  ipoints_3dmovie=0
  do ispec=1,NSPEC_CRUST_MANTLE
    do k=1,NGLLZ,iNIT
      do j=1,NGLLY,iNIT
        do i=1,NGLLX,iNIT
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
    if(mask_ibool(iglob)) then  !this element is in the region
      ispecele  = ispecele+1
      do k=1,NGLLZ-1,iNIT
        do j=1,NGLLY-1,iNIT
          do i=1,NGLLX-1,iNIT
            ! if(mask_3dmovie(i,j,k,ispec)) then
            iglob1 = ibool_crust_mantle(i,j,k,ispec)
            iglob2 = ibool_crust_mantle(i+iNIT,j,k,ispec)
            iglob3 = ibool_crust_mantle(i+iNIT,j+iNIT,k,ispec)
            iglob4 = ibool_crust_mantle(i,j+iNIT,k,ispec)
            iglob5 = ibool_crust_mantle(i,j,k+iNIT,ispec)
            iglob6 = ibool_crust_mantle(i+iNIT,j,k+iNIT,ispec)
            iglob7 = ibool_crust_mantle(i+iNIT,j+iNIT,k+iNIT,ispec)
            iglob8 = ibool_crust_mantle(i,j+iNIT,k+iNIT,ispec)
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

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_movie_volume_strains(myrank,npoints_3dmovie, &
                                        LOCAL_PATH,MOVIE_VOLUME_TYPE,MOVIE_COARSE, &
                    it,muvstore_crust_mantle_3dmovie,mask_3dmovie,nu_3dmovie,&
                    hprime_xx,hprime_yy,hprime_zz,ibool_crust_mantle,&
                    xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle,&
                    etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle,&
                    gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle,displ_crust_mantle)

! outputs strains: MOVIE_VOLUME_TYPE == 1 / 2 / 3

  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  ! input
  integer :: myrank,npoints_3dmovie,MOVIE_VOLUME_TYPE,it

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: eps_trace_over_3_loc_crust_mantle
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: &
        xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle,&
        etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
        gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_3DMOVIE) :: &
    muvstore_crust_mantle_3dmovie

  logical, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: mask_3dmovie

  logical :: MOVIE_COARSE
  real(kind=CUSTOM_REAL), dimension(3,3,npoints_3dmovie) :: nu_3dmovie
  character(len=150) LOCAL_PATH,outputname

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: displ_crust_mantle

  ! variables
  character(len=150) prname
  integer :: ipoints_3dmovie,i,j,k,ispec,iNIT,ier
  real(kind=CUSTOM_REAL) :: muv_3dmovie
  real(kind=CUSTOM_REAL),dimension(3,3) :: eps_loc,eps_loc_new
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: store_val3d_NN,store_val3d_EE,store_val3d_ZZ,&
                                                     store_val3d_NE,store_val3d_NZ,store_val3d_EZ

  character(len=1) movie_prefix

  ! check
  if(NDIM /= 3) call exit_MPI(myrank, 'write_movie_volume requires NDIM = 3')

  ! allocates arrays
  allocate(store_val3d_NN(npoints_3dmovie), &
          store_val3d_EE(npoints_3dmovie), &
          store_val3d_ZZ(npoints_3dmovie), &
          store_val3d_NE(npoints_3dmovie), &
          store_val3d_NZ(npoints_3dmovie), &
          store_val3d_EZ(npoints_3dmovie), &
          stat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error allocating store_val3d_ .. arrays')

  if(MOVIE_VOLUME_TYPE == 1) then
    movie_prefix='E' ! strain
  else if(MOVIE_VOLUME_TYPE == 2) then
    movie_prefix='S' ! time integral of strain
  else if(MOVIE_VOLUME_TYPE == 3) then
    movie_prefix='P' ! potency, or integral of strain x \mu
  endif
  if(MOVIE_COARSE) then
   iNIT = NGLLX-1
  else
   iNIT = 1
  endif

  write(prname,"('proc',i6.6)") myrank

  ipoints_3dmovie=0
  do ispec=1,NSPEC_CRUST_MANTLE
   call compute_element_strain_undo_att_noDev(ispec,nglob_crust_mantle,NSPEC_CRUST_MANTLE,&
                                              displ_crust_mantle,hprime_xx,hprime_yy,hprime_zz,ibool_crust_mantle,&
                                              xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle,&
                                              etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle,&
                                              gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle,&
                                              epsilondev_loc_crust_mantle,eps_trace_over_3_loc_crust_mantle)

   do k=1,NGLLZ,iNIT
    do j=1,NGLLY,iNIT
     do i=1,NGLLX,iNIT
      if(mask_3dmovie(i,j,k,ispec)) then
       ipoints_3dmovie=ipoints_3dmovie+1
       muv_3dmovie=muvstore_crust_mantle_3dmovie(i,j,k,ispec)

       eps_loc(1,1)=eps_trace_over_3_loc_crust_mantle(i,j,k) + epsilondev_loc_crust_mantle(1,i,j,k)
       eps_loc(2,2)=eps_trace_over_3_loc_crust_mantle(i,j,k) + epsilondev_loc_crust_mantle(2,i,j,k)
       eps_loc(3,3)=eps_trace_over_3_loc_crust_mantle(i,j,k)- &
                 epsilondev_loc_crust_mantle(1,i,j,k) - epsilondev_loc_crust_mantle(2,i,j,k)
       eps_loc(1,2)=epsilondev_loc_crust_mantle(3,i,j,k)
       eps_loc(1,3)=epsilondev_loc_crust_mantle(4,i,j,k)
       eps_loc(2,3)=epsilondev_loc_crust_mantle(5,i,j,k)

       eps_loc(2,1)=eps_loc(1,2)
       eps_loc(3,1)=eps_loc(1,3)
       eps_loc(3,2)=eps_loc(2,3)

       ! rotate eps_loc to spherical coordinates
       eps_loc_new(:,:) = matmul(matmul(nu_3dmovie(:,:,ipoints_3dmovie),eps_loc(:,:)), &
                                  transpose(nu_3dmovie(:,:,ipoints_3dmovie)))
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

  deallocate(store_val3d_NN,store_val3d_EE,store_val3d_ZZ, &
            store_val3d_NE,store_val3d_NZ,store_val3d_EZ)

  end subroutine write_movie_volume_strains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_movie_volume_vector(myrank,it,npoints_3dmovie,LOCAL_PATH,MOVIE_VOLUME_TYPE, &
                                      MOVIE_COARSE,ibool_crust_mantle,vector_crust_mantle, &
                                      scalingval,mask_3dmovie,nu_3dmovie)

! outputs displacement/velocity: MOVIE_VOLUME_TYPE == 5 / 6

  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  ! input
  integer :: myrank,it
  integer :: npoints_3dmovie
  integer :: MOVIE_VOLUME_TYPE
  logical :: MOVIE_COARSE

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  real(kind=CUSTOM_REAL), dimension(3,NGLOB_CRUST_MANTLE) :: vector_crust_mantle

  logical, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: mask_3dmovie

  double precision :: scalingval
  real(kind=CUSTOM_REAL), dimension(3,3,npoints_3dmovie) :: nu_3dmovie

  character(len=150) LOCAL_PATH

  real(kind=CUSTOM_REAL) :: scalingval_to_use

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(3) :: vector_local,vector_local_new
  integer :: ipoints_3dmovie,i,j,k,ispec,iNIT,iglob,ier
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: store_val3d_N,store_val3d_E,store_val3d_Z

  character(len=150) outputname
  character(len=2) movie_prefix

  ! check
  if(NDIM /= 3) call exit_MPI(myrank,'write_movie_volume requires NDIM = 3')

  ! allocates arrays
  allocate(store_val3d_N(npoints_3dmovie), &
          store_val3d_E(npoints_3dmovie), &
          store_val3d_Z(npoints_3dmovie), &
          stat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error allocating store_val3d_N,.. movie arrays')

  if(MOVIE_VOLUME_TYPE == 5) then
    movie_prefix='DI' ! displacement
  else if(MOVIE_VOLUME_TYPE == 6) then
    movie_prefix='VE' ! velocity
  endif

  if(MOVIE_COARSE) then
   iNIT = NGLLX-1
  else
   iNIT = 1
  endif

  if(CUSTOM_REAL == SIZE_REAL) then
    scalingval_to_use = sngl(scalingval)
  else
    scalingval_to_use = scalingval
  endif

  ipoints_3dmovie = 0

  do ispec=1,NSPEC_CRUST_MANTLE
   do k=1,NGLLZ,iNIT
    do j=1,NGLLY,iNIT
     do i=1,NGLLX,iNIT
      if(mask_3dmovie(i,j,k,ispec)) then
       ipoints_3dmovie=ipoints_3dmovie+1
       iglob = ibool_crust_mantle(i,j,k,ispec)
       vector_local(:) = vector_crust_mantle(:,iglob) * scalingval_to_use

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

  deallocate(store_val3d_N,store_val3d_E,store_val3d_Z)

  end subroutine write_movie_volume_vector

!
!-------------------------------------------------------------------------------------------------
!

 subroutine write_movie_volume_divcurl(myrank,it,eps_trace_over_3_crust_mantle,&
                                      div_displ_outer_core, &
                                      accel_outer_core,kappavstore_outer_core,rhostore_outer_core,ibool_outer_core, &
                                      eps_trace_over_3_inner_core, &
                                      epsilondev_crust_mantle,epsilondev_inner_core, &
                                      LOCAL_PATH, &
                                      displ_crust_mantle,displ_inner_core,displ_outer_core, &
                                      veloc_crust_mantle,veloc_inner_core,veloc_outer_core, &
                                      accel_crust_mantle,accel_inner_core, &
                                      ibool_crust_mantle,ibool_inner_core)

! outputs divergence and curl: MOVIE_VOLUME_TYPE == 4

  implicit none
  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer :: myrank,it
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: eps_trace_over_3_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ADJOINT) :: div_displ_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: accel_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: &
        rhostore_outer_core,kappavstore_outer_core
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: ibool_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_STRAIN_ONLY) :: eps_trace_over_3_inner_core
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: epsilondev_crust_mantle
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_STR_OR_ATT) :: epsilondev_inner_core

  character(len=150) LOCAL_PATH

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: displ_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: veloc_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: accel_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: veloc_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: displ_outer_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE) :: displ_inner_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE) :: veloc_inner_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE) :: accel_inner_core

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: ibool_inner_core

  ! local parameters
  real(kind=CUSTOM_REAL) :: rhol,kappal
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: div_s_outer_core
  integer :: ispec,iglob,i,j,k,ier
  character(len=150) outputname
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: tmp_data
  real(kind=CUSTOM_REAL) :: scale_displ,scale_veloc,scale_accel

  ! output parameters
  logical,parameter :: MOVIE_OUTPUT_DIV = .true.          ! divergence
  logical,parameter :: MOVIE_OUTPUT_CURL = .false.        ! curl
  logical,parameter :: MOVIE_OUTPUT_CURLNORM = .true.     ! Frobenius norm of curl
  logical,parameter :: MOVIE_OUTPUT_DISPLNORM = .false.   ! norm of displacement
  logical,parameter :: MOVIE_OUTPUT_VELOCNORM = .false.   ! norm of velocity
  logical,parameter :: MOVIE_OUTPUT_ACCELNORM = .false.   ! norm of acceleration

  ! outputs divergence
  if( MOVIE_OUTPUT_DIV ) then
    ! crust_mantle region
    ! these binary arrays can be converted into mesh format using the utilitiy ./bin/xcombine_vol_data
    ! old name format:     write(outputname,"('proc',i6.6,'_crust_mantle_div_displ_it',i6.6,'.bin')") myrank,it
    write(outputname,"('proc',i6.6,'_reg1_div_displ_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'error opening file '//trim(outputname))
    write(27) eps_trace_over_3_crust_mantle
    close(27)

    ! outer core
    if (NSPEC_OUTER_CORE_ADJOINT /= 1 ) then
      write(outputname,"('proc',i6.6,'_reg2_div_displ_it',i6.6,'.bin')") myrank,it
      open(unit=27,file=trim(LOCAL_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
      if( ier /= 0 ) call exit_MPI(myrank,'error opening file '//trim(outputname))
      write(27)  ONE_THIRD * div_displ_outer_core
      close(27)
    else
      ! we use div s = - p / kappa = rhostore_outer_core * accel_outer_core / kappavstore_outer_core
      allocate(div_s_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE))
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

      ! old name format: write(outputname,"('proc',i6.6,'_outer_core_div_displ_it',i6.6,'.bin')") myrank,it
      write(outputname,"('proc',i6.6,'_reg2_div_displ_it',i6.6,'.bin')") myrank,it
      open(unit=27,file=trim(LOCAL_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
      if( ier /= 0 ) call exit_MPI(myrank,'error opening file '//trim(outputname))
      write(27)  div_s_outer_core
      close(27)

      deallocate(div_s_outer_core)
    endif

    ! inner core
    ! old name format: write(outputname,"('proc',i6.6,'_inner_core_div_displ_proc_it',i6.6,'.bin')") myrank,it
    write(outputname,"('proc',i6.6,'_reg3_div_displ_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'error opening file '//trim(outputname))
    write(27) eps_trace_over_3_inner_core
    close(27)

  endif

  ! outputs epsilondev
  if( MOVIE_OUTPUT_CURL ) then
    ! epsilondev
    ! these binary files must be handled by the user, no further utilities available for this format
    ! crust mantle
    write(outputname,"('proc',i6.6,'_crust_mantle_epsdev_displ_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//trim(outputname),status='unknown',form='unformatted')
    write(27) epsilondev_crust_mantle
    close(27)
    ! inner core
    write(outputname,"('proc',i6.6,'inner_core_epsdev_displ_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//trim(outputname),status='unknown',form='unformatted')
    write(27) epsilondev_inner_core
    close(27)
  endif

  ! outputs norm of epsilondev
  if( MOVIE_OUTPUT_CURLNORM ) then
    ! these binary arrays can be converted into mesh format using the utilitiy ./bin/xcombine_vol_data
    ! crust_mantle region
    write(outputname,"('proc',i6.6,'_reg1_epsdev_displ_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'error opening file '//trim(outputname))
    ! frobenius norm
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE))
    do ispec = 1, NSPEC_CRUST_MANTLE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            tmp_data(i,j,k,ispec) = sqrt( epsilondev_crust_mantle(1,i,j,k,ispec)**2 &
                                          + epsilondev_crust_mantle(2,i,j,k,ispec)**2 &
                                          + epsilondev_crust_mantle(3,i,j,k,ispec)**2 &
                                          + epsilondev_crust_mantle(4,i,j,k,ispec)**2 &
                                          + epsilondev_crust_mantle(5,i,j,k,ispec)**2)
          enddo
        enddo
      enddo
    enddo
    write(27) tmp_data
    close(27)
    deallocate(tmp_data)

    ! alternative: e.g. first component only
    !write(27) epsilondev_crust_mantle(1,:,:,:,:)
    !close(27)

    ! inner core
    write(outputname,"('proc',i6.6,'_reg3_epsdev_displ_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'error opening file '//trim(outputname))
    ! frobenius norm
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE))
    do ispec = 1, NSPEC_INNER_CORE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            tmp_data(i,j,k,ispec) = sqrt( epsilondev_inner_core(1,i,j,k,ispec)**2 &
                                          + epsilondev_inner_core(2,i,j,k,ispec)**2 &
                                          + epsilondev_inner_core(3,i,j,k,ispec)**2 &
                                          + epsilondev_inner_core(4,i,j,k,ispec)**2 &
                                          + epsilondev_inner_core(5,i,j,k,ispec)**2)
          enddo
        enddo
      enddo
    enddo
    write(27) tmp_data
    close(27)
    deallocate(tmp_data)

    ! alternative: e.g. first component only
    !write(27) epsilondev_inner_core(1,:,:,:,:)
    !close(27)
  endif

  ! dimensionalized scalings
  scale_displ = R_EARTH
  scale_veloc = scale_displ * sqrt(PI*GRAV*RHOAV)
  scale_accel = scale_veloc * dsqrt(PI*GRAV*RHOAV)

  ! outputs norm of displacement
  if( MOVIE_OUTPUT_DISPLNORM ) then
    ! crust mantle
    ! these binary arrays can be converted into mesh format using the utilitiy ./bin/xcombine_vol_data
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE))
    do ispec = 1, NSPEC_CRUST_MANTLE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_crust_mantle(i,j,k,ispec)
            ! norm
            tmp_data(i,j,k,ispec) = scale_displ * sqrt( displ_crust_mantle(1,iglob)**2 &
                                          + displ_crust_mantle(2,iglob)**2 &
                                          + displ_crust_mantle(3,iglob)**2 )
          enddo
        enddo
      enddo
    enddo
    write(outputname,"('proc',i6.6,'_reg1_displ_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'error opening file '//trim(outputname))
    write(27) tmp_data
    close(27)
    deallocate(tmp_data)

    ! outer core
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE))
    do ispec = 1, NSPEC_OUTER_CORE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_outer_core(i,j,k,ispec)
            ! norm
            ! note: disp_outer_core is potential, this just outputs the potential,
            !          not the actual displacement u = grad(rho * Chi) / rho
            tmp_data(i,j,k,ispec) = abs(displ_outer_core(iglob))
          enddo
        enddo
      enddo
    enddo
    write(outputname,"('proc',i6.6,'_reg2_displ_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'error opening file '//trim(outputname))
    write(27) tmp_data
    close(27)
    deallocate(tmp_data)

    ! inner core
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE))
    do ispec = 1, NSPEC_INNER_CORE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_inner_core(i,j,k,ispec)
            ! norm
            tmp_data(i,j,k,ispec) = scale_displ * sqrt( displ_inner_core(1,iglob)**2 &
                                          + displ_inner_core(2,iglob)**2 &
                                          + displ_inner_core(3,iglob)**2 )
          enddo
        enddo
      enddo
    enddo
    write(outputname,"('proc',i6.6,'_reg3_displ_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'error opening file '//trim(outputname))
    write(27) tmp_data
    close(27)
    deallocate(tmp_data)
  endif


  ! outputs norm of velocity
  if( MOVIE_OUTPUT_VELOCNORM ) then
    ! crust mantle
    ! these binary arrays can be converted into mesh format using the utilitiy ./bin/xcombine_vol_data
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE))
    do ispec = 1, NSPEC_CRUST_MANTLE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_crust_mantle(i,j,k,ispec)
            ! norm of velocity
            tmp_data(i,j,k,ispec) = scale_veloc * sqrt( veloc_crust_mantle(1,iglob)**2 &
                                          + veloc_crust_mantle(2,iglob)**2 &
                                          + veloc_crust_mantle(3,iglob)**2 )
          enddo
        enddo
      enddo
    enddo
    write(outputname,"('proc',i6.6,'_reg1_veloc_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'error opening file '//trim(outputname))
    write(27) tmp_data
    close(27)
    deallocate(tmp_data)

    ! outer core
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE))
    do ispec = 1, NSPEC_OUTER_CORE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_outer_core(i,j,k,ispec)
            ! norm of velocity
            ! note: this outputs only the first time derivative of the potential,
            !          not the actual velocity v = grad(Chi_dot)
            tmp_data(i,j,k,ispec) = abs(veloc_outer_core(iglob))
          enddo
        enddo
      enddo
    enddo
    write(outputname,"('proc',i6.6,'_reg2_veloc_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'error opening file '//trim(outputname))
    write(27) tmp_data
    close(27)
    deallocate(tmp_data)

    ! inner core
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE))
    do ispec = 1, NSPEC_INNER_CORE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_inner_core(i,j,k,ispec)
            ! norm of velocity
            tmp_data(i,j,k,ispec) = scale_veloc * sqrt( veloc_inner_core(1,iglob)**2 &
                                          + veloc_inner_core(2,iglob)**2 &
                                          + veloc_inner_core(3,iglob)**2 )
          enddo
        enddo
      enddo
    enddo
    write(outputname,"('proc',i6.6,'_reg3_veloc_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'error opening file '//trim(outputname))
    write(27) tmp_data
    close(27)
    deallocate(tmp_data)
  endif

  ! outputs norm of acceleration
  if( MOVIE_OUTPUT_ACCELNORM ) then
    ! acceleration
    ! these binary arrays can be converted into mesh format using the utilitiy ./bin/xcombine_vol_data
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE))
    do ispec = 1, NSPEC_CRUST_MANTLE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_crust_mantle(i,j,k,ispec)
            ! norm
            tmp_data(i,j,k,ispec) = scale_accel * sqrt( accel_crust_mantle(1,iglob)**2 &
                                          + accel_crust_mantle(2,iglob)**2 &
                                          + accel_crust_mantle(3,iglob)**2 )
          enddo
        enddo
      enddo
    enddo
    write(outputname,"('proc',i6.6,'_reg1_accel_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'error opening file '//trim(outputname))
    write(27) tmp_data
    close(27)
    deallocate(tmp_data)

    ! outer core acceleration
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE))
    do ispec = 1, NSPEC_OUTER_CORE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_outer_core(i,j,k,ispec)
            ! norm
            ! note: this outputs only the second time derivative of the potential,
            !          not the actual acceleration or pressure p = - rho * Chi_dot_dot
            tmp_data(i,j,k,ispec) = abs(accel_outer_core(iglob))
          enddo
        enddo
      enddo
    enddo
    write(outputname,"('proc',i6.6,'_reg2_accel_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'error opening file '//trim(outputname))
    write(27) tmp_data
    close(27)
    deallocate(tmp_data)

    ! inner core
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE))
    do ispec = 1, NSPEC_INNER_CORE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_inner_core(i,j,k,ispec)
            ! norm of velocity
            tmp_data(i,j,k,ispec) = scale_accel * sqrt( accel_inner_core(1,iglob)**2 &
                                          + accel_inner_core(2,iglob)**2 &
                                          + accel_inner_core(3,iglob)**2 )
          enddo
        enddo
      enddo
    enddo
    write(outputname,"('proc',i6.6,'_reg3_accel_it',i6.6,'.bin')") myrank,it
    open(unit=27,file=trim(LOCAL_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'error opening file '//trim(outputname))
    write(27) tmp_data
    close(27)
    deallocate(tmp_data)
  endif

  end subroutine write_movie_volume_divcurl


!-------------------------------------------------------------------------------------------------
! external mesh routine for saving vtk files for custom_real values on global points

  subroutine write_VTK_data_cr(idoubling,nspec,nglob, &
                              xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                              glob_data,prname_file)

! outputs single file for each process

  implicit none

  include "constants.h"

  integer :: nspec,nglob

  integer, dimension(nspec):: idoubling

  ! global coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy

  ! global data values array
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob) :: glob_data

  ! file name
  character(len=256) prname_file

  ! local parameters
  integer :: ispec,i
  real(kind=CUSTOM_REAL) :: rval,thetaval,phival,xval,yval,zval

  ! write source and receiver VTK files for Paraview
  open(IOVTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown')
  write(IOVTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOVTK,'(a)') 'material model VTK file'
  write(IOVTK,'(a)') 'ASCII'
  write(IOVTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOVTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do i=1,nglob

    !x,y,z store have been converted to r theta phi already, need to revert back for xyz output
    rval = xstore_dummy(i)
    thetaval = ystore_dummy(i)
    phival = zstore_dummy(i)
    call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)

    !write(IOVTK,'(3e18.6)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)
    write(IOVTK,'(3e18.6)') xval,yval,zval
  enddo
  write(IOVTK,*) ! blank line

  ! defines cell on coarse corner points
  ! note: indices for vtk start at 0
  write(IOVTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
  do ispec=1,nspec

    ! specific to inner core elements
    ! exclude fictitious elements in central cube
    if(idoubling(ispec) /= IFLAG_IN_FICTITIOUS_CUBE) then
      ! valid cell
      write(IOVTK,'(9i12)') 8,ibool(1,1,1,ispec)-1, &
                          ibool(NGLLX,1,1,ispec)-1, &
                          ibool(NGLLX,NGLLY,1,ispec)-1, &
                          ibool(1,NGLLY,1,ispec)-1, &
                          ibool(1,1,NGLLZ,ispec)-1, &
                          ibool(NGLLX,1,NGLLZ,ispec)-1, &
                          ibool(NGLLX,NGLLY,NGLLZ,ispec)-1, &
                          ibool(1,NGLLY,NGLLZ,ispec)-1
    else
      ! fictitious elements in central cube
      ! maps cell onto a randomly chosen point
      write(IOVTK,'(9i12)') 8,ibool(1,1,1,1)-1, &
                            ibool(1,1,1,1)-1, &
                            ibool(1,1,1,1)-1, &
                            ibool(1,1,1,1)-1, &
                            ibool(1,1,1,1)-1, &
                            ibool(1,1,1,1)-1, &
                            ibool(1,1,1,1)-1, &
                            ibool(1,1,1,1)-1
    endif

  enddo
  write(IOVTK,*) ! blank line

  ! type: hexahedrons
  write(IOVTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOVTK,*) (12,ispec=1,nspec)
  write(IOVTK,*) ! blank line

  ! x components
  write(IOVTK,'(a,i12)') "POINT_DATA ",nglob
  write(IOVTK,'(a)') "SCALARS x_comp float"
  write(IOVTK,'(a)') "LOOKUP_TABLE default"
  do i = 1,nglob
      write(IOVTK,*) glob_data(1,i)
  enddo
  ! y components
  write(IOVTK,'(a)') "SCALARS y_comp float"
  write(IOVTK,'(a)') "LOOKUP_TABLE default"
  do i = 1,nglob
      write(IOVTK,*) glob_data(2,i)
  enddo
  ! z components
  write(IOVTK,'(a)') "SCALARS z_comp float"
  write(IOVTK,'(a)') "LOOKUP_TABLE default"
  do i = 1,nglob
      write(IOVTK,*) glob_data(3,i)
  enddo
  ! norm
  write(IOVTK,'(a)') "SCALARS norm float"
  write(IOVTK,'(a)') "LOOKUP_TABLE default"
  do i = 1,nglob
      write(IOVTK,*) sqrt( glob_data(1,i)*glob_data(1,i) &
                        + glob_data(2,i)*glob_data(2,i) &
                        + glob_data(3,i)*glob_data(3,i))
  enddo
  write(IOVTK,*) ! blank line

  close(IOVTK)


  end subroutine write_VTK_data_cr

!-------------------------------------------------------------------------------------------------

! external mesh routine for saving vtk files for custom_real values on global points

  subroutine write_VTK_data_cr_all(myrank,idoubling, &
                              nspec,nglob, &
                              xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                              glob_data,prname_file)

! outputs single file for all processes

  use mpi

  implicit none

  include "precision.h"
  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer :: myrank,nspec,nglob

  integer, dimension(nspec):: idoubling

  ! global coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy

  ! global data values array
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob) :: glob_data

  ! file name
  character(len=256) prname_file

  ! local parameters
  integer :: ispec,i,iproc,ier
  real(kind=CUSTOM_REAL) :: rval,thetaval,phival,xval,yval,zval

  real(kind=CUSTOM_REAL), dimension(:,:),allocatable :: &
      store_val_x_all,store_val_y_all,store_val_z_all, &
      store_val_ux_all,store_val_uy_all,store_val_uz_all
  integer, dimension(:,:,:,:,:),allocatable :: ibool_all
  integer, dimension(:,:),allocatable :: idoubling_all

  ! master collect arrays
  if( myrank == 0 ) then
    allocate(store_val_x_all(nglob,0:NPROCTOT_VAL-1), &
          store_val_y_all(nglob,0:NPROCTOT_VAL-1), &
          store_val_z_all(nglob,0:NPROCTOT_VAL-1), &
          store_val_ux_all(nglob,0:NPROCTOT_VAL-1), &
          store_val_uy_all(nglob,0:NPROCTOT_VAL-1), &
          store_val_uz_all(nglob,0:NPROCTOT_VAL-1), &
          idoubling_all(nspec,0:NPROCTOT_VAL-1), &
          ibool_all(NGLLX,NGLLY,NGLLZ,nspec,0:NPROCTOT_VAL-1),stat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error allocating stores')
  endif

  ! gather info on master proc
  call MPI_GATHER(xstore_dummy,nglob,CUSTOM_MPI_TYPE,store_val_x_all,nglob,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(ystore_dummy,nglob,CUSTOM_MPI_TYPE,store_val_y_all,nglob,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(zstore_dummy,nglob,CUSTOM_MPI_TYPE,store_val_z_all,nglob,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)

  call MPI_GATHER(glob_data(1,:),nglob,CUSTOM_MPI_TYPE,store_val_ux_all,nglob,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(glob_data(2,:),nglob,CUSTOM_MPI_TYPE,store_val_uy_all,nglob,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(glob_data(3,:),nglob,CUSTOM_MPI_TYPE,store_val_uz_all,nglob,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)

  call MPI_GATHER(ibool,NGLLX*NGLLY*NGLLZ*nspec,MPI_INTEGER,ibool_all, &
                  NGLLX*NGLLY*NGLLZ*nspec,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(idoubling,nspec,MPI_INTEGER,idoubling_all,nspec,MPI_INTEGER,0,MPI_COMM_WORLD,ier)


  if( myrank == 0 ) then

    ! write source and receiver VTK files for Paraview
    open(IOVTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown')
    write(IOVTK,'(a)') '# vtk DataFile Version 3.1'
    write(IOVTK,'(a)') 'material model VTK file'
    write(IOVTK,'(a)') 'ASCII'
    write(IOVTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
    write(IOVTK, '(a,i12,a)') 'POINTS ', nglob*NPROCTOT_VAL, ' float'
    do iproc=0, NPROCTOT_VAL-1
      do i=1,nglob

        !x,y,z store have been converted to r theta phi already, need to revert back for xyz output
        rval = store_val_x_all(i,iproc)
        thetaval = store_val_y_all(i,iproc)
        phival = store_val_z_all(i,iproc)
        call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)

        !write(IOVTK,'(3e18.6)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)
        write(IOVTK,'(3e18.6)') xval,yval,zval
      enddo
    enddo
    write(IOVTK,*) ! blank line

    ! defines cell on coarse corner points
    ! note: indices for vtk start at 0
    write(IOVTK,'(a,i12,i12)') "CELLS ",nspec*NPROCTOT_VAL,nspec*NPROCTOT_VAL*9
    do iproc=0, NPROCTOT_VAL-1
      do ispec=1,nspec

        ! note: central cube elements are only shared and used in CHUNK_AB and CHUNK_AB_ANTIPODE
        !          all other chunks ignore those elements

        ! specific to inner core elements
        ! exclude fictitious elements in central cube
        if(idoubling_all(ispec,iproc) /= IFLAG_IN_FICTITIOUS_CUBE) then
          ! valid cell
          ! cell corner ids
          write(IOVTK,'(9i12)') 8,ibool_all(1,1,1,ispec,iproc)-1+iproc*nglob, &
                            ibool_all(NGLLX,1,1,ispec,iproc)-1+iproc*nglob, &
                            ibool_all(NGLLX,NGLLY,1,ispec,iproc)-1+iproc*nglob, &
                            ibool_all(1,NGLLY,1,ispec,iproc)-1+iproc*nglob, &
                            ibool_all(1,1,NGLLZ,ispec,iproc)-1+iproc*nglob, &
                            ibool_all(NGLLX,1,NGLLZ,ispec,iproc)-1+iproc*nglob, &
                            ibool_all(NGLLX,NGLLY,NGLLZ,ispec,iproc)-1+iproc*nglob, &
                            ibool_all(1,NGLLY,NGLLZ,ispec,iproc)-1+iproc*nglob
        else
          ! fictitious elements in central cube
          ! maps cell onto a randomly chosen point
          write(IOVTK,'(9i12)') 8,ibool_all(1,1,1,1,iproc)-1, &
                            ibool_all(1,1,1,1,iproc)-1, &
                            ibool_all(1,1,1,1,iproc)-1, &
                            ibool_all(1,1,1,1,iproc)-1, &
                            ibool_all(1,1,1,1,iproc)-1, &
                            ibool_all(1,1,1,1,iproc)-1, &
                            ibool_all(1,1,1,1,iproc)-1, &
                            ibool_all(1,1,1,1,iproc)-1
        endif

      enddo
    enddo
    write(IOVTK,*) ! blank line

    ! type: hexahedrons
    write(IOVTK,'(a,i12)') "CELL_TYPES ",nspec*NPROCTOT_VAL
    write(IOVTK,*) (12,ispec=1,nspec*NPROCTOT_VAL)
    write(IOVTK,*) ! blank line

    ! x components
    write(IOVTK,'(a,i12)') "POINT_DATA ",nglob*NPROCTOT_VAL
    write(IOVTK,'(a)') "SCALARS x_comp float"
    write(IOVTK,'(a)') "LOOKUP_TABLE default"
    do iproc=0, NPROCTOT_VAL-1
      do i = 1,nglob
        write(IOVTK,*) store_val_ux_all(i,iproc)
      enddo
    enddo
    ! y components
    write(IOVTK,'(a)') "SCALARS y_comp float"
    write(IOVTK,'(a)') "LOOKUP_TABLE default"
    do iproc=0, NPROCTOT_VAL-1
      do i = 1,nglob
        write(IOVTK,*) store_val_uy_all(i,iproc)
      enddo
    enddo
    ! z components
    write(IOVTK,'(a)') "SCALARS z_comp float"
    write(IOVTK,'(a)') "LOOKUP_TABLE default"
    do iproc=0, NPROCTOT_VAL-1
      do i = 1,nglob
        write(IOVTK,*) store_val_uz_all(i,iproc)
      enddo
    enddo
    ! norm
    write(IOVTK,'(a)') "SCALARS norm float"
    write(IOVTK,'(a)') "LOOKUP_TABLE default"
    do iproc=0, NPROCTOT_VAL-1
      do i = 1,nglob
        write(IOVTK,*) sqrt( store_val_ux_all(i,iproc)**2 &
                          + store_val_uy_all(i,iproc)**2 &
                          + store_val_uz_all(i,iproc)**2 )
      enddo
    enddo
    write(IOVTK,*) ! blank line

    close(IOVTK)

  endif

  if( myrank == 0 ) then
    deallocate(store_val_x_all,store_val_y_all,store_val_z_all, &
            store_val_ux_all,store_val_uy_all,store_val_uz_all, &
            ibool_all)
  endif

  end subroutine write_VTK_data_cr_all
