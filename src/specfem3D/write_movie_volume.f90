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

! this subroutine counts the number of points and elements within the movie volume
! in this processor slice, and returns arrays that keep track of them, both in global and local indexing schemes

  subroutine movie_volume_count_points(num_ibool_3dmovie,mask_ibool_3dmovie)

  use specfem_par, only: CUSTOM_REAL,MAX_STRING_LEN,NGLLX,NGLLY,NGLLZ,IOUT,NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE_3DMOVIE, &
    LOCAL_TMP_PATH,MOVIE_COARSE,MOVIE_TOP,MOVIE_BOTTOM,MOVIE_EAST,MOVIE_WEST,MOVIE_NORTH,MOVIE_SOUTH, &
    myrank

  use specfem_par_crustmantle, only: ibool_crust_mantle,rstore_crust_mantle

  use specfem_par_movie, only: mask_3dmovie,nspecel_3dmovie,npoints_3dmovie

  implicit none

  integer, dimension(NGLOB_CRUST_MANTLE_3DMOVIE),intent(out) :: num_ibool_3dmovie
  logical, dimension(NGLOB_CRUST_MANTLE_3DMOVIE),intent(out) :: mask_ibool_3dmovie

  ! local parameters
  integer :: ipoints_3dmovie,ispecel_3dmovie
  integer :: ispec,iglob,iglob_center,i,j,k,iNIT,ier
  real(kind=CUSTOM_REAL) :: rval,thetaval,phival
  character(len=MAX_STRING_LEN) :: prname_movie

  ! initializes movie points
  num_ibool_3dmovie(:) = -99
  mask_ibool_3dmovie(:) = .false.
  mask_3dmovie(:,:,:,:) = .false.

  ! stepping
  if (MOVIE_COARSE) then
    iNIT = NGLLX-1
  else
    iNIT = 1
  endif

  !find and count points within given region for storing movie
  ipoints_3dmovie = 0
  ispecel_3dmovie = 0
  do ispec = 1,NSPEC_CRUST_MANTLE
    ! store if center of element is in the given region
    iglob_center = ibool_crust_mantle((NGLLX+1)/2,(NGLLY+1)/2,(NGLLZ+1)/2,ispec)

    rval     = rstore_crust_mantle(1,iglob_center)
    thetaval = rstore_crust_mantle(2,iglob_center)
    phival   = rstore_crust_mantle(3,iglob_center)

    ! we already changed xyz back to rthetaphi
    if ( (rval < MOVIE_TOP .and. rval > MOVIE_BOTTOM) .and. &
         (thetaval > MOVIE_NORTH .and. thetaval < MOVIE_SOUTH) .and. &
         ( (phival < MOVIE_EAST .and. phival > MOVIE_WEST) .or. &
         ( (MOVIE_EAST < MOVIE_WEST) .and. (phival > MOVIE_EAST .or. phival < MOVIE_WEST) ) )) then

      ! adds as movie element
      ispecel_3dmovie = ispecel_3dmovie + 1

      ! sets flag at center for later checks
      if (MOVIE_COARSE) mask_ibool_3dmovie(iglob_center) = .true.

      do k = 1,NGLLZ,iNIT
        do j = 1,NGLLY,iNIT
          do i = 1,NGLLX,iNIT
            iglob = ibool_crust_mantle(i,j,k,ispec)

            if (.not. mask_ibool_3dmovie(iglob)) then
              ! stores point
              ipoints_3dmovie = ipoints_3dmovie + 1

              mask_ibool_3dmovie(iglob) = .true.
              mask_3dmovie(i,j,k,ispec) = .true.
              num_ibool_3dmovie(iglob) = ipoints_3dmovie
            endif

          enddo !i
        enddo !j
      enddo !k

    endif !in region
  enddo !ispec
  npoints_3dmovie = ipoints_3dmovie
  nspecel_3dmovie = ispecel_3dmovie

  ! create name of database
  write(prname_movie,'(a,i6.6,a)') trim(LOCAL_TMP_PATH)//'/'//'proc',myrank,'_'
  open(unit=IOUT,file=trim(prname_movie)//'movie3D_info.txt', &
        status='unknown',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening file movie3D_info.txt')

  write(IOUT,*) npoints_3dmovie, nspecel_3dmovie
  close(IOUT)

  end subroutine movie_volume_count_points

!
!-------------------------------------------------------------------------------------------------
!


! put the list of parameters back here to avoid a warning / error from the gfortran compiler
! about undefined behavior when aggressive loop vectorization is used by the compiler

  subroutine write_movie_volume_mesh(nu_3dmovie,num_ibool_3dmovie,mask_3dmovie,mask_ibool_3dmovie, &
                                     muvstore_crust_mantle_3dmovie,npoints_3dmovie)

! writes meshfiles to merge with solver snapshots for 3D volume movies.  Also computes and outputs
! the rotation matrix nu_3dmovie required to transfer to a geographic coordinate system

  use specfem_par
  use specfem_par_crustmantle, only: ibool_crust_mantle,rstore_crust_mantle

  implicit none

  integer,intent(in) :: npoints_3dmovie
  integer, dimension(NGLOB_CRUST_MANTLE_3DMOVIE),intent(in) :: num_ibool_3dmovie

  real(kind=CUSTOM_REAL), dimension(3,3,npoints_3dmovie),intent(inout) :: nu_3dmovie

  logical, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_3DMOVIE),intent(in) :: mask_3dmovie
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_3DMOVIE),intent(in) :: muvstore_crust_mantle_3dmovie
  logical, dimension(NGLOB_CRUST_MANTLE_3DMOVIE),intent(in) :: mask_ibool_3dmovie

  ! local parameters
  integer :: ipoints_3dmovie,ispecele,ispec,i,j,k,iNIT,ier
  integer :: iglob,iglob_center
  integer :: iglob1,iglob2,iglob3,iglob4,iglob5,iglob6,iglob7,iglob8
  integer :: n1,n2,n3,n4,n5,n6,n7,n8
  real(kind=CUSTOM_REAL) :: rval,thetaval,phival,xval,yval,zval,st,ct,sp,cp
  real(kind=CUSTOM_REAL), dimension(npoints_3dmovie) :: store_val3D_x,store_val3D_y, store_val3D_z
  real(kind=CUSTOM_REAL), dimension(npoints_3dmovie) :: store_val3D_mu
  character(len=MAX_STRING_LEN) :: prname_movie

  ! outputs values into ascii file for each process/slice
  logical,parameter :: DEBUG_OUTPUT_ASCII = .false.

  ! safety check
  if (NDIM /= 3) stop 'movie volume output requires NDIM = 3'

  ! stepping
  if (MOVIE_COARSE) then
    iNIT = NGLLX-1
  else
    iNIT = 1
  endif

  ! loops over all elements
  ipoints_3dmovie = 0
  do ispec = 1,NSPEC_CRUST_MANTLE

    ! checks center of element for movie flag
    iglob_center = ibool_crust_mantle((NGLLX+1)/2,(NGLLY+1)/2,(NGLLZ+1)/2,ispec)

    ! checks if movie element
    if (mask_ibool_3dmovie(iglob_center)) then

      ! stores element coordinates
      do k = 1,NGLLZ,iNIT
        do j = 1,NGLLY,iNIT
          do i = 1,NGLLX,iNIT
            ! only store points once
            if (mask_3dmovie(i,j,k,ispec)) then
              ! point increment
              ipoints_3dmovie = ipoints_3dmovie + 1

              ! gets point position
              iglob = ibool_crust_mantle(i,j,k,ispec)

              rval = rstore_crust_mantle(1,iglob)
              thetaval = rstore_crust_mantle(2,iglob)
              phival = rstore_crust_mantle(3,iglob)

              !x,y,z store have been converted to r theta phi already, need to revert back for xyz output
              call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)

              store_val3D_x(ipoints_3dmovie) = xval
              store_val3D_y(ipoints_3dmovie) = yval
              store_val3D_z(ipoints_3dmovie) = zval
              store_val3D_mu(ipoints_3dmovie) = muvstore_crust_mantle_3dmovie(i,j,k,ispec)

              st = sin(thetaval)
              ct = cos(thetaval)
              sp = sin(phival)
              cp = cos(phival)

              nu_3dmovie(1,1,ipoints_3dmovie) = -ct*cp
              nu_3dmovie(1,2,ipoints_3dmovie) = -ct*sp
              nu_3dmovie(1,3,ipoints_3dmovie) = st
              nu_3dmovie(2,1,ipoints_3dmovie) = -sp
              nu_3dmovie(2,2,ipoints_3dmovie) = cp
              nu_3dmovie(2,3,ipoints_3dmovie) = 0.d0
              nu_3dmovie(3,1,ipoints_3dmovie) = st*cp
              nu_3dmovie(3,2,ipoints_3dmovie) = st*sp
              nu_3dmovie(3,3,ipoints_3dmovie) = ct
            endif !mask_3dmovie
          enddo  !i
        enddo  !j
      enddo  !k
    endif

  enddo !ispec

  write(prname_movie,'(a,i6.6,a)') trim(LOCAL_TMP_PATH)//'/'//'proc',myrank,'_'

  ! stores point locations
  open(unit=IOUT,file=trim(prname_movie)//'movie3D_x.bin',status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening file movie3D_x.bin')
  if (npoints_3dmovie > 0) then
    write(IOUT) store_val3D_x(1:npoints_3dmovie)
  endif
  close(IOUT)

  open(unit=IOUT,file=trim(prname_movie)//'movie3D_y.bin',status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening file movie3D_y.bin')
  if (npoints_3dmovie > 0) then
    write(IOUT) store_val3D_y(1:npoints_3dmovie)
   endif
  close(IOUT)

  open(unit=IOUT,file=trim(prname_movie)//'movie3D_z.bin',status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening file movie3D_z.bin')
  if (npoints_3dmovie > 0) then
    write(IOUT) store_val3D_z(1:npoints_3dmovie)
  endif
  close(IOUT)

  if (DEBUG_OUTPUT_ASCII) then
    open(unit=IOUT,file=trim(prname_movie)//'ascii_output.txt',status='unknown',iostat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error opening file ascii_output.txt')
    if (npoints_3dmovie > 0) then
      do i = 1,npoints_3dmovie
        write(IOUT,*) store_val3D_x(i),store_val3D_y(i),store_val3D_z(i),store_val3D_mu(i)
      enddo
    endif
    close(IOUT)
  endif

  ! store element infos
  open(unit=IOUT,file=trim(prname_movie)//'movie3D_elements.bin',status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening file movie3D_elements.bin')
 !  open(unit=57,file=trim(prname_movie)//'movie3D_elements.txt',status='unknown')

  ispecele = 0
  do ispec = 1,NSPEC_CRUST_MANTLE

    ! checks center of element for movie flag
    iglob_center = ibool_crust_mantle((NGLLX+1)/2,(NGLLY+1)/2,(NGLLZ+1)/2,ispec)

    ! checks if movie element
    if (mask_ibool_3dmovie(iglob_center)) then

      ! this element is in the movie region
      ispecele = ispecele+1

      do k = 1,NGLLZ-1,iNIT
        do j = 1,NGLLY-1,iNIT
          do i = 1,NGLLX-1,iNIT
            ! defines corners of a vtk element
            iglob1 = ibool_crust_mantle(i,j,k,ispec)
            iglob2 = ibool_crust_mantle(i+iNIT,j,k,ispec)
            iglob3 = ibool_crust_mantle(i+iNIT,j+iNIT,k,ispec)
            iglob4 = ibool_crust_mantle(i,j+iNIT,k,ispec)
            iglob5 = ibool_crust_mantle(i,j,k+iNIT,ispec)
            iglob6 = ibool_crust_mantle(i+iNIT,j,k+iNIT,ispec)
            iglob7 = ibool_crust_mantle(i+iNIT,j+iNIT,k+iNIT,ispec)
            iglob8 = ibool_crust_mantle(i,j+iNIT,k+iNIT,ispec)

            ! vtk indexing starts at 0 -> adds minus 1
            n1 = num_ibool_3dmovie(iglob1)-1
            n2 = num_ibool_3dmovie(iglob2)-1
            n3 = num_ibool_3dmovie(iglob3)-1
            n4 = num_ibool_3dmovie(iglob4)-1
            n5 = num_ibool_3dmovie(iglob5)-1
            n6 = num_ibool_3dmovie(iglob6)-1
            n7 = num_ibool_3dmovie(iglob7)-1
            n8 = num_ibool_3dmovie(iglob8)-1

            ! checks indices
            if (n1 < 0 .or. n2 < 0 .or. n3 < 0 .or. n4 < 0 .or. n5 < 0 .or. n6 < 0 .or. n7 < 0 .or. n8 < 0) then
              print *,'Error: movie element ',ispec,ispecele,'has invalid node index:',n1,n2,n3,n4,n5,n6,n7,n8
              call exit_mpi(myrank,'Error invalid movie element node index')
            endif

            write(IOUT) n1,n2,n3,n4,n5,n6,n7,n8

            ! text output
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

  subroutine movie_volume_integrate_strain(deltat,vnspec, &
                                           eps_trace_over_3, &
                                           epsilondev_xx,epsilondev_yy, &
                                           epsilondev_xy,epsilondev_xz, &
                                           epsilondev_yz, &
                                           Ieps_trace_over_3, &
                                           Iepsilondev_xx,Iepsilondev_yy, &
                                           Iepsilondev_xy,Iepsilondev_xz, &
                                           Iepsilondev_yz)

  use constants_solver

  implicit none

  real(kind=CUSTOM_REAL),intent(in) :: deltat

  integer,intent(in) :: vnspec
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,vnspec),intent(in) :: &
    eps_trace_over_3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,vnspec),intent(in) :: &
    epsilondev_xx,epsilondev_yy,epsilondev_xy, &
    epsilondev_xz,epsilondev_yz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,vnspec),intent(inout) :: &
    Ieps_trace_over_3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,vnspec),intent(inout) :: &
    Iepsilondev_xx,Iepsilondev_yy,Iepsilondev_xy, &
    Iepsilondev_xz,Iepsilondev_yz

  ! updates integral values
  Iepsilondev_xx(:,:,:,:) = Iepsilondev_xx(:,:,:,:) + deltat*epsilondev_xx(:,:,:,:)
  Iepsilondev_yy(:,:,:,:) = Iepsilondev_yy(:,:,:,:) + deltat*epsilondev_yy(:,:,:,:)
  Iepsilondev_xy(:,:,:,:) = Iepsilondev_xy(:,:,:,:) + deltat*epsilondev_xy(:,:,:,:)
  Iepsilondev_xz(:,:,:,:) = Iepsilondev_xz(:,:,:,:) + deltat*epsilondev_xz(:,:,:,:)
  Iepsilondev_yz(:,:,:,:) = Iepsilondev_yz(:,:,:,:) + deltat*epsilondev_yz(:,:,:,:)

  Ieps_trace_over_3(:,:,:,:) = Ieps_trace_over_3(:,:,:,:) + deltat*eps_trace_over_3(:,:,:,:)

  end subroutine movie_volume_integrate_strain

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_movie_volume_strains(vnspec_eps_cm, &
                                        eps_trace_over_3_crust_mantle, &
                                        vnspec_cm, &
                                        epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
                                        epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle)

! outputs strains: MOVIE_VOLUME_TYPE == 1 / 2 / 3

  use constants_solver

  use shared_parameters, only: LOCAL_TMP_PATH,MOVIE_VOLUME_TYPE,MOVIE_COARSE

  use specfem_par, only: it

  use specfem_par_movie, only: npoints_3dmovie,muvstore_crust_mantle_3dmovie,mask_3dmovie,nu_3dmovie

  implicit none

  ! input
  integer,intent(in) :: vnspec_eps_cm
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,vnspec_eps_cm),intent(in) :: eps_trace_over_3_crust_mantle

  integer,intent(in) :: vnspec_cm
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,vnspec_cm),intent(in) :: &
    epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
    epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle

  ! variables
  real(kind=CUSTOM_REAL) :: muv_3dmovie
  real(kind=CUSTOM_REAL),dimension(3,3) :: eps_loc,eps_loc_new
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: store_val3d_NN,store_val3d_EE,store_val3d_ZZ, &
                                                     store_val3d_NE,store_val3d_NZ,store_val3d_EZ
  integer :: ipoints_3dmovie,i,j,k,ispec,iNIT,ier
  character(len=MAX_STRING_LEN) :: outputname
  character(len=1) :: movie_prefix

  ! check
  if (NDIM /= 3) call exit_MPI(myrank, 'write_movie_volume_strains() requires NDIM = 3')
  if (vnspec_cm /= NSPEC_CRUST_MANTLE) call exit_MPI(myrank,'Invalid vnspec_cm value for write_movie_volume_strains() routine')

  ! allocates arrays
  allocate(store_val3d_NN(npoints_3dmovie), &
           store_val3d_EE(npoints_3dmovie), &
           store_val3d_ZZ(npoints_3dmovie), &
           store_val3d_NE(npoints_3dmovie), &
           store_val3d_NZ(npoints_3dmovie), &
           store_val3d_EZ(npoints_3dmovie), &
           stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating store_val3d_ .. arrays')

  if (MOVIE_VOLUME_TYPE == 1) then
    movie_prefix='E' ! strain
  else if (MOVIE_VOLUME_TYPE == 2) then
    movie_prefix='S' ! time integral of strain
  else if (MOVIE_VOLUME_TYPE == 3) then
    movie_prefix='P' ! potency, or integral of strain x \mu
  endif

  ! stepping
  if (MOVIE_COARSE) then
   iNIT = NGLLX-1
  else
   iNIT = 1
  endif

  ipoints_3dmovie = 0
  do ispec = 1,NSPEC_CRUST_MANTLE
    do k = 1,NGLLZ,iNIT
      do j = 1,NGLLY,iNIT
        do i = 1,NGLLX,iNIT
          if (mask_3dmovie(i,j,k,ispec)) then
            ipoints_3dmovie = ipoints_3dmovie + 1
            muv_3dmovie = muvstore_crust_mantle_3dmovie(i,j,k,ispec)

            eps_loc(1,1) = eps_trace_over_3_crust_mantle(i,j,k,ispec) + epsilondev_xx_crust_mantle(i,j,k,ispec)
            eps_loc(2,2) = eps_trace_over_3_crust_mantle(i,j,k,ispec) + epsilondev_yy_crust_mantle(i,j,k,ispec)
            eps_loc(3,3) = eps_trace_over_3_crust_mantle(i,j,k,ispec) &
                           - epsilondev_xx_crust_mantle(i,j,k,ispec) &
                           - epsilondev_yy_crust_mantle(i,j,k,ispec)

            eps_loc(1,2) = epsilondev_xy_crust_mantle(i,j,k,ispec)
            eps_loc(1,3) = epsilondev_xz_crust_mantle(i,j,k,ispec)
            eps_loc(2,3) = epsilondev_yz_crust_mantle(i,j,k,ispec)

            eps_loc(2,1) = eps_loc(1,2)
            eps_loc(3,1) = eps_loc(1,3)
            eps_loc(3,2) = eps_loc(2,3)

            ! rotate eps_loc to spherical coordinates
            eps_loc_new(:,:) = matmul(matmul(nu_3dmovie(:,:,ipoints_3dmovie),eps_loc(:,:)), &
                                      transpose(nu_3dmovie(:,:,ipoints_3dmovie)))
            if (MOVIE_VOLUME_TYPE == 3) eps_loc_new(:,:) = eps_loc(:,:)*muv_3dmovie

            store_val3d_NN(ipoints_3dmovie) = eps_loc_new(1,1)
            store_val3d_EE(ipoints_3dmovie) = eps_loc_new(2,2)
            store_val3d_ZZ(ipoints_3dmovie) = eps_loc_new(3,3)
            store_val3d_NE(ipoints_3dmovie) = eps_loc_new(1,2)
            store_val3d_NZ(ipoints_3dmovie) = eps_loc_new(1,3)
            store_val3d_EZ(ipoints_3dmovie) = eps_loc_new(2,3)
          endif
        enddo
      enddo
    enddo
  enddo
  if (ipoints_3dmovie /= npoints_3dmovie) stop 'did not find the right number of points for 3D movie'

  write(outputname,"('proc',i6.6,'_movie3D_',a,'NN',i6.6,'.bin')") myrank,movie_prefix,it
  open(unit=IOUT,file=trim(LOCAL_TMP_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening file '//trim(outputname))
  write(IOUT) store_val3d_NN(1:npoints_3dmovie)
  close(IOUT)

  write(outputname,"('proc',i6.6,'_movie3D_',a,'EE',i6.6,'.bin')") myrank,movie_prefix,it
  open(unit=IOUT,file=trim(LOCAL_TMP_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening file '//trim(outputname))
  write(IOUT) store_val3d_EE(1:npoints_3dmovie)
  close(IOUT)

  write(outputname,"('proc',i6.6,'_movie3D_',a,'ZZ',i6.6,'.bin')") myrank,movie_prefix,it
  open(unit=IOUT,file=trim(LOCAL_TMP_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening file '//trim(outputname))
  write(IOUT) store_val3d_ZZ(1:npoints_3dmovie)
  close(IOUT)

  write(outputname,"('proc',i6.6,'_movie3D_',a,'NE',i6.6,'.bin')") myrank,movie_prefix,it
  open(unit=IOUT,file=trim(LOCAL_TMP_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening file '//trim(outputname))
  write(IOUT) store_val3d_NE(1:npoints_3dmovie)
  close(IOUT)

  write(outputname,"('proc',i6.6,'_movie3D_',a,'NZ',i6.6,'.bin')") myrank,movie_prefix,it
  open(unit=IOUT,file=trim(LOCAL_TMP_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening file '//trim(outputname))
  write(IOUT) store_val3d_NZ(1:npoints_3dmovie)
  close(IOUT)

  write(outputname,"('proc',i6.6,'_movie3D_',a,'EZ',i6.6,'.bin')") myrank,movie_prefix,it
  open(unit=IOUT,file=trim(LOCAL_TMP_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening file '//trim(outputname))
  write(IOUT) store_val3d_EZ(1:npoints_3dmovie)
  close(IOUT)

  deallocate(store_val3d_NN,store_val3d_EE,store_val3d_ZZ, &
             store_val3d_NE,store_val3d_NZ,store_val3d_EZ)

  end subroutine write_movie_volume_strains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_movie_volume_vector(npoints_3dmovie, &
                                       ibool_crust_mantle,vector_crust_mantle, &
                                       scalingval,mask_3dmovie,nu_3dmovie)

! outputs displacement/velocity: MOVIE_VOLUME_TYPE == 5 / 6

  use constants_solver

  use shared_parameters, only: LOCAL_TMP_PATH,MOVIE_VOLUME_TYPE,MOVIE_COARSE

  use specfem_par, only: it

  implicit none

  ! input
  integer :: npoints_3dmovie

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle

  ! displacement or velocity array
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: vector_crust_mantle

  logical, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_3DMOVIE) :: mask_3dmovie

  double precision :: scalingval
  real(kind=CUSTOM_REAL), dimension(NDIM,NDIM,npoints_3dmovie) :: nu_3dmovie

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(NDIM) :: vector_local,vector_local_new
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: store_val3d_N,store_val3d_E,store_val3d_Z

  integer :: ipoints_3dmovie,i,j,k,ispec,iNIT,iglob,ier
  character(len=MAX_STRING_LEN) :: outputname
  character(len=2) :: movie_prefix

  ! check
  if (NDIM /= 3) call exit_MPI(myrank,'write_movie_volume requires NDIM = 3')

  ! allocates arrays
  allocate(store_val3d_N(npoints_3dmovie), &
           store_val3d_E(npoints_3dmovie), &
           store_val3d_Z(npoints_3dmovie), &
           stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating store_val3d_N,.. movie arrays')

  if (MOVIE_VOLUME_TYPE == 5) then
    movie_prefix='DI' ! displacement
  else if (MOVIE_VOLUME_TYPE == 6) then
    movie_prefix='VE' ! velocity
  endif

  if (MOVIE_COARSE) then
   iNIT = NGLLX-1
  else
   iNIT = 1
  endif

  ipoints_3dmovie = 0

  ! stores field in crust/mantle region
  do ispec = 1,NSPEC_CRUST_MANTLE
    do k = 1,NGLLZ,iNIT
      do j = 1,NGLLY,iNIT
        do i = 1,NGLLX,iNIT
          if (mask_3dmovie(i,j,k,ispec)) then
            ipoints_3dmovie = ipoints_3dmovie + 1
            iglob = ibool_crust_mantle(i,j,k,ispec)

            ! dimensionalizes field by scaling
            vector_local(:) = vector_crust_mantle(:,iglob)*real(scalingval,kind=CUSTOM_REAL)

            ! rotate eps_loc to spherical coordinates
            vector_local_new(:) = matmul(nu_3dmovie(:,:,ipoints_3dmovie), vector_local(:))

            ! stores field
            store_val3d_N(ipoints_3dmovie) = vector_local_new(1)
            store_val3d_E(ipoints_3dmovie) = vector_local_new(2)
            store_val3d_Z(ipoints_3dmovie) = vector_local_new(3)
          endif
        enddo
      enddo
   enddo
  enddo
  close(IOUT)

  ! checks number of processed points
  if (ipoints_3dmovie /= npoints_3dmovie) stop 'did not find the right number of points for 3D movie'

  ! file output
  write(outputname,"('proc',i6.6,'_movie3D_',a,'N',i6.6,'.bin')") myrank,movie_prefix,it
  open(unit=IOUT,file=trim(LOCAL_TMP_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening file '//trim(outputname))
  write(IOUT) store_val3d_N(1:npoints_3dmovie)
  close(IOUT)

  write(outputname,"('proc',i6.6,'_movie3D_',a,'E',i6.6,'.bin')") myrank,movie_prefix,it
  open(unit=IOUT,file=trim(LOCAL_TMP_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening file '//trim(outputname))
  write(IOUT) store_val3d_E(1:npoints_3dmovie)
  close(IOUT)

  write(outputname,"('proc',i6.6,'_movie3D_',a,'Z',i6.6,'.bin')") myrank,movie_prefix,it
  open(unit=IOUT,file=trim(LOCAL_TMP_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening file '//trim(outputname))
  write(IOUT) store_val3d_Z(1:npoints_3dmovie)
  close(IOUT)

  deallocate(store_val3d_N,store_val3d_E,store_val3d_Z)

  end subroutine write_movie_volume_vector

!
!-------------------------------------------------------------------------------------------------
!

 subroutine write_movie_volume_divcurl(vnspec_eps_cm,eps_trace_over_3_crust_mantle, &
                                       div_displ_outer_core, &
                                       accel_outer_core,kappavstore_outer_core,rhostore_outer_core,ibool_outer_core, &
                                       vnspec_eps_ic,eps_trace_over_3_inner_core, &
                                       vnspec_cm,epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
                                       epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle, &
                                       vnspec_ic,epsilondev_xx_inner_core,epsilondev_yy_inner_core,epsilondev_xy_inner_core, &
                                       epsilondev_xz_inner_core,epsilondev_yz_inner_core)

! outputs divergence and curl: MOVIE_VOLUME_TYPE == 4

  use constants_solver

  use shared_parameters, only: LOCAL_TMP_PATH

  use specfem_par, only: it

  implicit none

  integer,intent(in) :: vnspec_eps_cm
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,vnspec_eps_cm) :: eps_trace_over_3_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_3DMOVIE) :: div_displ_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: accel_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: rhostore_outer_core,kappavstore_outer_core
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: ibool_outer_core

  integer,intent(in) :: vnspec_eps_ic
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,vnspec_eps_ic) :: eps_trace_over_3_inner_core

  integer,intent(in) :: vnspec_cm
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,vnspec_cm) :: &
    epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
    epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle

  integer,intent(in) :: vnspec_ic
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,vnspec_ic) :: &
    epsilondev_xx_inner_core,epsilondev_yy_inner_core,epsilondev_xy_inner_core, &
    epsilondev_xz_inner_core,epsilondev_yz_inner_core

  ! local parameters
  real(kind=CUSTOM_REAL) :: rhol,kappal
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: div_s_outer_core
  integer :: ispec,iglob,i,j,k,ier
  character(len=MAX_STRING_LEN) outputname
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: tmp_data

  ! output parameters
  logical,parameter :: MOVIE_OUTPUT_DIV = .true.          ! divergence
  logical,parameter :: MOVIE_OUTPUT_CURL = .true.         ! curl
  logical,parameter :: MOVIE_OUTPUT_CURLNORM = .true.     ! Frobenius norm of curl

  ! checks
  if (vnspec_cm /= NSPEC_CRUST_MANTLE) call exit_MPI(myrank,'Invalid vnspec_cm value for write_movie_volume_divcurl() routine')

  ! outputs divergence
  if (MOVIE_OUTPUT_DIV) then
    ! crust_mantle region
    ! these binary arrays can be converted into mesh format using the utility ./bin/xcombine_vol_data
    ! old name format:     write(outputname,"('proc',i6.6,'_crust_mantle_div_displ_it',i6.6,'.bin')") myrank,it
    write(outputname,"('proc',i6.6,'_reg1_div_displ_it',i6.6,'.bin')") myrank,it
    open(unit=IOUT,file=trim(LOCAL_TMP_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening file '//trim(outputname))
    write(IOUT) eps_trace_over_3_crust_mantle
    close(IOUT)

    ! outer core
    if (NSPEC_OUTER_CORE_3DMOVIE > 1) then
      write(outputname,"('proc',i6.6,'_reg2_div_displ_it',i6.6,'.bin')") myrank,it
      open(unit=IOUT,file=trim(LOCAL_TMP_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error opening file '//trim(outputname))
      write(IOUT)  ONE_THIRD * div_displ_outer_core
      close(IOUT)
    else
      ! we use div s = - p / kappa = rhostore_outer_core * accel_outer_core / kappavstore_outer_core
      allocate(div_s_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array div_s_outer_core')

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
      open(unit=IOUT,file=trim(LOCAL_TMP_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error opening file '//trim(outputname))
      write(IOUT)  div_s_outer_core
      close(IOUT)

      deallocate(div_s_outer_core)
    endif

    ! inner core
    ! old name format: write(outputname,"('proc',i6.6,'_inner_core_div_displ_proc_it',i6.6,'.bin')") myrank,it
    write(outputname,"('proc',i6.6,'_reg3_div_displ_it',i6.6,'.bin')") myrank,it
    open(unit=IOUT,file=trim(LOCAL_TMP_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening file '//trim(outputname))
    write(IOUT) eps_trace_over_3_inner_core
    close(IOUT)

  endif

  ! outputs epsilondev
  if (MOVIE_OUTPUT_CURL) then
    ! epsilondev
    ! these binary files must be handled by the user, no further utilities available for this format
    ! crust mantle
    write(outputname,"('proc',i6.6,'_crust_mantle_epsdev_displ_it',i6.6,'.bin')") myrank,it
    open(unit=IOUT,file=trim(LOCAL_TMP_PATH)//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error opening file '//trim(outputname))
    write(IOUT) epsilondev_xx_crust_mantle
    write(IOUT) epsilondev_yy_crust_mantle
    write(IOUT) epsilondev_xy_crust_mantle
    write(IOUT) epsilondev_xz_crust_mantle
    write(IOUT) epsilondev_yz_crust_mantle
    close(IOUT)
    ! inner core
    write(outputname,"('proc',i6.6,'inner_core_epsdev_displ_it',i6.6,'.bin')") myrank,it
    open(unit=IOUT,file=trim(LOCAL_TMP_PATH)//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error opening file '//trim(outputname))
    write(IOUT) epsilondev_xx_inner_core
    write(IOUT) epsilondev_yy_inner_core
    write(IOUT) epsilondev_xy_inner_core
    write(IOUT) epsilondev_xz_inner_core
    write(IOUT) epsilondev_yz_inner_core
    close(IOUT)
  endif

  ! outputs norm of epsilondev
  if (MOVIE_OUTPUT_CURLNORM) then
    ! these binary arrays can be converted into mesh format using the utility ./bin/xcombine_vol_data
    ! crust_mantle region
    write(outputname,"('proc',i6.6,'_reg1_epsdev_displ_it',i6.6,'.bin')") myrank,it
    open(unit=IOUT,file=trim(LOCAL_TMP_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening file '//trim(outputname))
    ! Frobenius norm
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')

    do ispec = 1, NSPEC_CRUST_MANTLE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            tmp_data(i,j,k,ispec) = sqrt( epsilondev_xx_crust_mantle(i,j,k,ispec)**2 &
                                          + epsilondev_yy_crust_mantle(i,j,k,ispec)**2 &
                                          + epsilondev_xy_crust_mantle(i,j,k,ispec)**2 &
                                          + epsilondev_xz_crust_mantle(i,j,k,ispec)**2 &
                                          + epsilondev_yz_crust_mantle(i,j,k,ispec)**2)
          enddo
        enddo
      enddo
    enddo
    write(IOUT) tmp_data
    close(IOUT)
    deallocate(tmp_data)

    ! inner core
    write(outputname,"('proc',i6.6,'_reg3_epsdev_displ_it',i6.6,'.bin')") myrank,it
    open(unit=IOUT,file=trim(LOCAL_TMP_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening file '//trim(outputname))
    ! Frobenius norm
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')

    do ispec = 1, NSPEC_INNER_CORE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            tmp_data(i,j,k,ispec) = sqrt( epsilondev_xx_inner_core(i,j,k,ispec)**2 &
                                          + epsilondev_yy_inner_core(i,j,k,ispec)**2 &
                                          + epsilondev_xy_inner_core(i,j,k,ispec)**2 &
                                          + epsilondev_xz_inner_core(i,j,k,ispec)**2 &
                                          + epsilondev_yz_inner_core(i,j,k,ispec)**2)
          enddo
        enddo
      enddo
    enddo
    write(IOUT) tmp_data
    close(IOUT)
    deallocate(tmp_data)
  endif

  end subroutine write_movie_volume_divcurl

!
!-------------------------------------------------------------------------------------------------
!

 subroutine write_movie_volume_displnorm(displ_crust_mantle,displ_inner_core,displ_outer_core, &
                                         ibool_crust_mantle,ibool_inner_core,ibool_outer_core)

! outputs norm of displacement: MOVIE_VOLUME_TYPE == 7

  use constants_solver

  use shared_parameters, only: LOCAL_TMP_PATH

  use specfem_par, only: it, scale_displ

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: displ_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE) :: displ_inner_core
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: displ_outer_core

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: ibool_inner_core
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: ibool_outer_core

  ! local parameters
  integer :: ispec,iglob,i,j,k,ier
  character(len=MAX_STRING_LEN) outputname
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: tmp_data

  logical,parameter :: OUTPUT_CRUST_MANTLE = .true.
  logical,parameter :: OUTPUT_OUTER_CORE = .true.
  logical,parameter :: OUTPUT_INNER_CORE = .true.

  ! outputs norm of displacement
  if (OUTPUT_CRUST_MANTLE) then
    ! crust mantle
    ! these binary arrays can be converted into mesh format using the utility ./bin/xcombine_vol_data
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')

    do ispec = 1, NSPEC_CRUST_MANTLE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_crust_mantle(i,j,k,ispec)
            ! norm
            tmp_data(i,j,k,ispec) = real(scale_displ,kind=CUSTOM_REAL) * sqrt( displ_crust_mantle(1,iglob)**2 &
                                          + displ_crust_mantle(2,iglob)**2 &
                                          + displ_crust_mantle(3,iglob)**2 )
          enddo
        enddo
      enddo
    enddo
    write(outputname,"('proc',i6.6,'_reg1_displ_it',i6.6,'.bin')") myrank,it
    open(unit=IOUT,file=trim(LOCAL_TMP_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening file '//trim(outputname))
    write(IOUT) tmp_data
    close(IOUT)
    deallocate(tmp_data)
  endif

  if (OUTPUT_OUTER_CORE) then
    ! outer core
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')

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
    open(unit=IOUT,file=trim(LOCAL_TMP_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening file '//trim(outputname))
    write(IOUT) tmp_data
    close(IOUT)
    deallocate(tmp_data)
  endif

  if (OUTPUT_INNER_CORE) then
    ! inner core
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')

    do ispec = 1, NSPEC_INNER_CORE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_inner_core(i,j,k,ispec)
            ! norm
            tmp_data(i,j,k,ispec) = real(scale_displ,kind=CUSTOM_REAL) * sqrt( displ_inner_core(1,iglob)**2 &
                                          + displ_inner_core(2,iglob)**2 &
                                          + displ_inner_core(3,iglob)**2 )
          enddo
        enddo
      enddo
    enddo
    write(outputname,"('proc',i6.6,'_reg3_displ_it',i6.6,'.bin')") myrank,it
    open(unit=IOUT,file=trim(LOCAL_TMP_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening file '//trim(outputname))
    write(IOUT) tmp_data
    close(IOUT)
    deallocate(tmp_data)
  endif

  end subroutine write_movie_volume_displnorm


!
!-------------------------------------------------------------------------------------------------
!

 subroutine write_movie_volume_velnorm(veloc_crust_mantle,veloc_inner_core,veloc_outer_core, &
                                       ibool_crust_mantle,ibool_inner_core,ibool_outer_core)

! outputs norm of velocity: MOVIE_VOLUME_TYPE == 8

  use constants_solver

  use shared_parameters, only: LOCAL_TMP_PATH

  use specfem_par, only: it, scale_veloc

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: veloc_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: veloc_outer_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE) :: veloc_inner_core

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: ibool_inner_core
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: ibool_outer_core

  ! local parameters
  integer :: ispec,iglob,i,j,k,ier
  character(len=MAX_STRING_LEN) :: outputname
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: tmp_data

  logical,parameter :: OUTPUT_CRUST_MANTLE = .true.
  logical,parameter :: OUTPUT_OUTER_CORE = .true.
  logical,parameter :: OUTPUT_INNER_CORE = .true.

  ! outputs norm of velocity
  if (OUTPUT_CRUST_MANTLE) then
    ! crust mantle
    ! these binary arrays can be converted into mesh format using the utility ./bin/xcombine_vol_data
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')

    do ispec = 1, NSPEC_CRUST_MANTLE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_crust_mantle(i,j,k,ispec)
            ! norm of velocity
            tmp_data(i,j,k,ispec) = real(scale_veloc,kind=CUSTOM_REAL) * sqrt( veloc_crust_mantle(1,iglob)**2 &
                                          + veloc_crust_mantle(2,iglob)**2 &
                                          + veloc_crust_mantle(3,iglob)**2 )
          enddo
        enddo
      enddo
    enddo
    write(outputname,"('proc',i6.6,'_reg1_veloc_it',i6.6,'.bin')") myrank,it
    open(unit=IOUT,file=trim(LOCAL_TMP_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening file '//trim(outputname))
    write(IOUT) tmp_data
    close(IOUT)
    deallocate(tmp_data)
  endif

  if (OUTPUT_OUTER_CORE) then
    ! outer core
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')

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
    open(unit=IOUT,file=trim(LOCAL_TMP_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening file '//trim(outputname))
    write(IOUT) tmp_data
    close(IOUT)
    deallocate(tmp_data)
  endif

  if (OUTPUT_INNER_CORE) then
    ! inner core
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')

    do ispec = 1, NSPEC_INNER_CORE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_inner_core(i,j,k,ispec)
            ! norm of velocity
            tmp_data(i,j,k,ispec) = real(scale_veloc,kind=CUSTOM_REAL) * sqrt( veloc_inner_core(1,iglob)**2 &
                                          + veloc_inner_core(2,iglob)**2 &
                                          + veloc_inner_core(3,iglob)**2 )
          enddo
        enddo
      enddo
    enddo
    write(outputname,"('proc',i6.6,'_reg3_veloc_it',i6.6,'.bin')") myrank,it
    open(unit=IOUT,file=trim(LOCAL_TMP_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening file '//trim(outputname))
    write(IOUT) tmp_data
    close(IOUT)
    deallocate(tmp_data)
  endif

  end subroutine write_movie_volume_velnorm

!
!-------------------------------------------------------------------------------------------------
!

 subroutine write_movie_volume_accelnorm(accel_crust_mantle,accel_inner_core,accel_outer_core, &
                                         ibool_crust_mantle,ibool_inner_core,ibool_outer_core)

! outputs norm of acceleration: MOVIE_VOLUME_TYPE == 1

  use constants_solver

  use shared_parameters, only: LOCAL_TMP_PATH

  use specfem_par, only: it, scale_t_inv,scale_veloc

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: accel_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE) :: accel_inner_core
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: accel_outer_core

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: ibool_inner_core
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: ibool_outer_core

  ! local parameters
  integer :: ispec,iglob,i,j,k,ier
  character(len=MAX_STRING_LEN) :: outputname
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: tmp_data
  real(kind=CUSTOM_REAL) :: scale_accel

  logical,parameter :: OUTPUT_CRUST_MANTLE = .true.
  logical,parameter :: OUTPUT_OUTER_CORE = .true.
  logical,parameter :: OUTPUT_INNER_CORE = .true.

  ! dimensionalized scaling
  scale_accel = real(scale_veloc * scale_t_inv,kind=CUSTOM_REAL)

  ! outputs norm of acceleration
  if (OUTPUT_CRUST_MANTLE) then
    ! acceleration
    ! these binary arrays can be converted into mesh format using the utility ./bin/xcombine_vol_data
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')

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
    open(unit=IOUT,file=trim(LOCAL_TMP_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening file '//trim(outputname))
    write(IOUT) tmp_data
    close(IOUT)
    deallocate(tmp_data)
  endif

  if (OUTPUT_OUTER_CORE) then
    ! outer core acceleration
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')

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
    open(unit=IOUT,file=trim(LOCAL_TMP_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening file '//trim(outputname))
    write(IOUT) tmp_data
    close(IOUT)
    deallocate(tmp_data)
  endif

  if (OUTPUT_INNER_CORE) then
    ! inner core
    allocate(tmp_data(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_data')

    do ispec = 1, NSPEC_INNER_CORE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_inner_core(i,j,k,ispec)
            ! norm of acceleration
            tmp_data(i,j,k,ispec) = scale_accel * sqrt( accel_inner_core(1,iglob)**2 &
                                          + accel_inner_core(2,iglob)**2 &
                                          + accel_inner_core(3,iglob)**2 )
          enddo
        enddo
      enddo
    enddo
    write(outputname,"('proc',i6.6,'_reg3_accel_it',i6.6,'.bin')") myrank,it
    open(unit=IOUT,file=trim(LOCAL_TMP_PATH)//'/'//trim(outputname),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening file '//trim(outputname))
    write(IOUT) tmp_data
    close(IOUT)
    deallocate(tmp_data)
  endif

  end subroutine write_movie_volume_accelnorm



