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


  subroutine write_VTK_data_points(nglob, &
                                  xstore_dummy,ystore_dummy,zstore_dummy, &
                                  points_globalindices,num_points_globalindices, &
                                  prname_file)

! external mesh routine for saving VTK files for points locations

  use constants,only: CUSTOM_REAL,MAX_STRING_LEN,IOUT_VTK

  implicit none

  integer,intent(in) :: nglob

  ! global coordinates
  real(kind=CUSTOM_REAL), dimension(nglob),intent(in) :: xstore_dummy,ystore_dummy,zstore_dummy

  ! GLL data values array
  integer,intent(in) :: num_points_globalindices
  integer, dimension(num_points_globalindices),intent(in) :: points_globalindices

  ! file name
  character(len=MAX_STRING_LEN),intent(in) :: prname_file

  integer :: i,iglob,ier

  ! write source and receiver VTK files for Paraview
  !debug
  !write(IMAIN,*) '  VTK file: '
  !write(IMAIN,*) '    ',prname_file(1:len_trim(prname_file))//'.vtk'

  open(IOUT_VTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown',iostat=ier)
  if (ier /= 0) stop 'Error opening VTK file'

  write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOUT_VTK,'(a)') 'material model VTK file'
  write(IOUT_VTK,'(a)') 'ASCII'
  write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOUT_VTK, '(a,i12,a)') 'POINTS ', num_points_globalindices, ' float'
  do i = 1,num_points_globalindices
    iglob = points_globalindices(i)
    if (iglob <= 0 .or. iglob > nglob) then
      print*,'Error: '//prname_file(1:len_trim(prname_file))//'.vtk'
      print*,'Error global index: ',iglob,i
      close(IOUT_VTK)
      stop 'Error VTK points file'
    endif

    write(IOUT_VTK,'(3e18.6)') xstore_dummy(iglob),ystore_dummy(iglob),zstore_dummy(iglob)
  enddo
  write(IOUT_VTK,*) ""

  close(IOUT_VTK)

  end subroutine write_VTK_data_points

!
!-------------------------------------------------------------------------------------------------
!


  subroutine write_VTK_glob_points(nglob, &
                                  xstore_dummy,ystore_dummy,zstore_dummy, &
                                  glob_values, &
                                  prname_file)

! external mesh routine for saving VTK files for points locations

  use constants,only: CUSTOM_REAL,MAX_STRING_LEN,IOUT_VTK

  implicit none

  integer,intent(in) :: nglob

  ! global coordinates
  real(kind=CUSTOM_REAL), dimension(nglob),intent(in) :: xstore_dummy,ystore_dummy,zstore_dummy

  ! GLL data values array
  real(kind=CUSTOM_REAL), dimension(nglob),intent(in) :: glob_values

  ! file name
  character(len=MAX_STRING_LEN),intent(in) :: prname_file

  ! local parameters
  integer :: iglob,ier

  ! write source and receiver VTK files for Paraview
  !debug
  !write(IMAIN,*) '  VTK file: '
  !write(IMAIN,*) '    ',prname_file(1:len_trim(prname_file))//'.vtk'

  open(IOUT_VTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown',iostat=ier)
  if (ier /= 0) stop 'Error opening VTK file'

  write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOUT_VTK,'(a)') 'material model VTK file'
  write(IOUT_VTK,'(a)') 'ASCII'
  write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOUT_VTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do iglob = 1,nglob
    write(IOUT_VTK,*) xstore_dummy(iglob),ystore_dummy(iglob),zstore_dummy(iglob)
  enddo
  write(IOUT_VTK,*) ""

  ! writes out GLL data (velocity) for each element point
  write(IOUT_VTK,'(a,i12)') "POINT_DATA ",nglob
  write(IOUT_VTK,'(a)') "SCALARS glob_data float"
  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
  do iglob = 1,nglob
    write(IOUT_VTK,*) glob_values(iglob)
  enddo
  write(IOUT_VTK,*) ""

  close(IOUT_VTK)

  end subroutine write_VTK_glob_points



!
!-------------------------------------------------------------------------------------------------
!


  subroutine write_VTK_data_elem_l(nspec,nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                        elem_flag,prname_file)

! routine for saving VTK file holding logical flag on each spectral element

  use constants,only: CUSTOM_REAL,MAX_STRING_LEN,NGLLX,NGLLY,NGLLZ,IOUT_VTK

  implicit none

  integer,intent(in) :: nspec,nglob

  ! global coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob),intent(in) :: xstore_dummy,ystore_dummy,zstore_dummy

  ! element flag array
  logical, dimension(nspec),intent(in) :: elem_flag

  ! file name
  character(len=MAX_STRING_LEN),intent(in) :: prname_file

  ! local parameters
  integer :: ispec,i,ier

  ! write source and receiver VTK files for Paraview
  !debug
  !write(IMAIN,*) '  VTK file: '
  !write(IMAIN,*) '    ',prname_file(1:len_trim(prname_file))//'.vtk'

  open(IOUT_VTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown',iostat=ier)
  if (ier /= 0) stop 'Error opening VTK file'

  write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOUT_VTK,'(a)') 'material model VTK file'
  write(IOUT_VTK,'(a)') 'ASCII'
  write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOUT_VTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do i = 1,nglob
    write(IOUT_VTK,'(3e18.6)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)
  enddo
  write(IOUT_VTK,*) ""

  ! note: indices for VTK start at 0
  write(IOUT_VTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
  do ispec = 1,nspec
    write(IOUT_VTK,'(9i12)') 8, &
          ibool(1,1,1,ispec)-1,ibool(NGLLX,1,1,ispec)-1,ibool(NGLLX,NGLLY,1,ispec)-1,ibool(1,NGLLY,1,ispec)-1,&
          ibool(1,1,NGLLZ,ispec)-1,ibool(NGLLX,1,NGLLZ,ispec)-1,ibool(NGLLX,NGLLY,NGLLZ,ispec)-1,ibool(1,NGLLY,NGLLZ,ispec)-1
  enddo
  write(IOUT_VTK,*) ""

  ! type: hexahedrons
  write(IOUT_VTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOUT_VTK,'(6i12)') (12,ispec = 1,nspec)
  write(IOUT_VTK,*) ""

  write(IOUT_VTK,'(a,i12)') "CELL_DATA ",nspec
  write(IOUT_VTK,'(a)') "SCALARS elem_flag integer"
  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
  do ispec = 1,nspec
    if (elem_flag(ispec) .eqv. .true.) then
      write(IOUT_VTK,*) 1
    else
      write(IOUT_VTK,*) 0
    endif
  enddo
  write(IOUT_VTK,*) ""
  close(IOUT_VTK)


  end subroutine write_VTK_data_elem_l


!
!-------------------------------------------------------------------------------------------------
!


  subroutine write_VTK_data_elem_i(nspec,nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                        elem_flag,prname_file)


! routine for saving VTK file holding integer value on each spectral element

  use constants,only: CUSTOM_REAL,MAX_STRING_LEN,NGLLX,NGLLY,NGLLZ,IOUT_VTK

  implicit none

  integer,intent(in) :: nspec,nglob

  ! global coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob),intent(in) :: xstore_dummy,ystore_dummy,zstore_dummy

  ! element flag array
  integer, dimension(nspec),intent(in) :: elem_flag

  ! file name
  character(len=MAX_STRING_LEN),intent(in) :: prname_file

  ! local parameters
  integer :: ispec,i,ier

  ! write source and receiver VTK files for Paraview
  !debug
  !write(IMAIN,*) '  VTK file: '
  !write(IMAIN,*) '    ',prname_file(1:len_trim(prname_file))//'.vtk'

  open(IOUT_VTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown',iostat=ier)
  if (ier /= 0) stop 'Error opening VTK file'

  write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOUT_VTK,'(a)') 'material model VTK file'
  write(IOUT_VTK,'(a)') 'ASCII'
  write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOUT_VTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do i = 1,nglob
    write(IOUT_VTK,'(3e18.6)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)
  enddo
  write(IOUT_VTK,*) ""

  ! note: indices for VTK start at 0
  write(IOUT_VTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
  do ispec = 1,nspec
    write(IOUT_VTK,'(9i12)') 8, &
          ibool(1,1,1,ispec)-1,ibool(NGLLX,1,1,ispec)-1,ibool(NGLLX,NGLLY,1,ispec)-1,ibool(1,NGLLY,1,ispec)-1,&
          ibool(1,1,NGLLZ,ispec)-1,ibool(NGLLX,1,NGLLZ,ispec)-1,ibool(NGLLX,NGLLY,NGLLZ,ispec)-1,ibool(1,NGLLY,NGLLZ,ispec)-1
  enddo
  write(IOUT_VTK,*) ""

  ! type: hexahedrons
  write(IOUT_VTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOUT_VTK,'(6i12)') (12,ispec = 1,nspec)
  write(IOUT_VTK,*) ""

  write(IOUT_VTK,'(a,i12)') "CELL_DATA ",nspec
  write(IOUT_VTK,'(a)') "SCALARS elem_val integer"
  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
  do ispec = 1,nspec
    write(IOUT_VTK,*) elem_flag(ispec)
  enddo
  write(IOUT_VTK,*) ""
  close(IOUT_VTK)

  end subroutine write_VTK_data_elem_i

!
!-------------------------------------------------------------------------------------------------
!

! external mesh routine for saving VTK files for custom_real values on global points

  subroutine write_VTK_data_cr(idoubling,nspec,nglob, &
                              xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                              glob_data,prname_file)

! outputs single file for each process

  use constants,only: CUSTOM_REAL,MAX_STRING_LEN,NDIM,NGLLX,NGLLY,NGLLZ,IOUT_VTK,IFLAG_IN_FICTITIOUS_CUBE

  implicit none

  integer,intent(in) :: nspec,nglob

  integer, dimension(nspec),intent(in) :: idoubling

  ! global coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob),intent(in) :: xstore_dummy,ystore_dummy,zstore_dummy

  ! global data values array
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(in) :: glob_data

  ! file name
  character(len=MAX_STRING_LEN),intent(in) :: prname_file

  ! local parameters
  integer :: ispec,i,ier
  real(kind=CUSTOM_REAL) :: rval,thetaval,phival,xval,yval,zval

  ! write source and receiver VTK files for Paraview
  open(IOUT_VTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown',iostat=ier)
  if (ier /= 0) stop 'Error opening VTK file'

  write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOUT_VTK,'(a)') 'material model VTK file'
  write(IOUT_VTK,'(a)') 'ASCII'
  write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOUT_VTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do i = 1,nglob

    !x,y,z store have been converted to r theta phi already, need to revert back for xyz output
    rval = xstore_dummy(i)
    thetaval = ystore_dummy(i)
    phival = zstore_dummy(i)
    call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)

    write(IOUT_VTK,'(3e18.6)') xval,yval,zval
  enddo
  write(IOUT_VTK,*) ""

  ! defines cell on coarse corner points
  ! note: indices for VTK start at 0
  write(IOUT_VTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
  do ispec = 1,nspec

    ! specific to inner core elements
    ! exclude fictitious elements in central cube
    if (idoubling(ispec) /= IFLAG_IN_FICTITIOUS_CUBE) then
      ! valid cell
      write(IOUT_VTK,'(9i12)') 8,ibool(1,1,1,ispec)-1, &
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
      write(IOUT_VTK,'(9i12)') 8,ibool(1,1,1,1)-1, &
                            ibool(1,1,1,1)-1, &
                            ibool(1,1,1,1)-1, &
                            ibool(1,1,1,1)-1, &
                            ibool(1,1,1,1)-1, &
                            ibool(1,1,1,1)-1, &
                            ibool(1,1,1,1)-1, &
                            ibool(1,1,1,1)-1
    endif

  enddo
  write(IOUT_VTK,*) ""

  ! type: hexahedrons
  write(IOUT_VTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOUT_VTK,'(6i12)') (12,ispec = 1,nspec)
  write(IOUT_VTK,*) ""

  ! x components
  write(IOUT_VTK,'(a,i12)') "POINT_DATA ",nglob
  write(IOUT_VTK,'(a)') "SCALARS x_comp float"
  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
  do i = 1,nglob
      write(IOUT_VTK,*) glob_data(1,i)
  enddo
  ! y components
  write(IOUT_VTK,'(a)') "SCALARS y_comp float"
  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
  do i = 1,nglob
      write(IOUT_VTK,*) glob_data(2,i)
  enddo
  ! z components
  write(IOUT_VTK,'(a)') "SCALARS z_comp float"
  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
  do i = 1,nglob
      write(IOUT_VTK,*) glob_data(3,i)
  enddo
  ! norm
  write(IOUT_VTK,'(a)') "SCALARS norm float"
  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
  do i = 1,nglob
      write(IOUT_VTK,*) sqrt( glob_data(1,i)*glob_data(1,i) &
                        + glob_data(2,i)*glob_data(2,i) &
                        + glob_data(3,i)*glob_data(3,i))
  enddo
  write(IOUT_VTK,*) ""

  close(IOUT_VTK)


  end subroutine write_VTK_data_cr

!
!-------------------------------------------------------------------------------------------------
!

! external mesh routine for saving VTK files for custom_real values on global points

  subroutine write_VTK_data_cr_all(myrank,NPROCTOT,idoubling, &
                              nspec,nglob, &
                              xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                              glob_data,prname_file)

! outputs single file for all processes

  use constants,only: CUSTOM_REAL,MAX_STRING_LEN,NDIM,NGLLX,NGLLY,NGLLZ,IOUT_VTK,IFLAG_IN_FICTITIOUS_CUBE

  implicit none

  integer,intent(in) :: myrank,NPROCTOT

  integer,intent(in) ::nspec,nglob

  integer, dimension(nspec),intent(in) :: idoubling

  ! global coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob),intent(in) :: xstore_dummy,ystore_dummy,zstore_dummy

  ! global data values array
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(in) :: glob_data

  ! file name
  character(len=MAX_STRING_LEN),intent(in) :: prname_file

  ! local parameters
  integer :: ispec,i,iproc,ier
  real(kind=CUSTOM_REAL) :: rval,thetaval,phival,xval,yval,zval

  real(kind=CUSTOM_REAL), dimension(:,:),allocatable :: &
      store_val_x_all,store_val_y_all,store_val_z_all, &
      store_val_ux_all,store_val_uy_all,store_val_uz_all
  integer, dimension(:,:,:,:,:),allocatable :: ibool_all
  integer, dimension(:,:),allocatable :: idoubling_all
  real(kind=CUSTOM_REAL), dimension(nglob) :: tmp

  ! master collect arrays
  if (myrank == 0) then
    allocate(store_val_x_all(nglob,0:NPROCTOT-1), &
            store_val_y_all(nglob,0:NPROCTOT-1), &
            store_val_z_all(nglob,0:NPROCTOT-1), &
            store_val_ux_all(nglob,0:NPROCTOT-1), &
            store_val_uy_all(nglob,0:NPROCTOT-1), &
            store_val_uz_all(nglob,0:NPROCTOT-1), &
            idoubling_all(nspec,0:NPROCTOT-1), &
            ibool_all(NGLLX,NGLLY,NGLLZ,nspec,0:NPROCTOT-1),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating stores')
  else
    ! dummy arrays
    allocate(store_val_x_all(1,1), &
            store_val_y_all(1,1), &
            store_val_z_all(1,1), &
            store_val_ux_all(1,1), &
            store_val_uy_all(1,1), &
            store_val_uz_all(1,1), &
            idoubling_all(1,1), &
            ibool_all(1,1,1,1,1),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating dummy stores')
  endif

  ! gather info on master proc
  call gather_all_cr(xstore_dummy,nglob,store_val_x_all,nglob,NPROCTOT)
  call gather_all_cr(ystore_dummy,nglob,store_val_y_all,nglob,NPROCTOT)
  call gather_all_cr(zstore_dummy,nglob,store_val_z_all,nglob,NPROCTOT)

  ! attention: these calls produce copies of the glob_data array
  ! x-component
  tmp(:) = glob_data(1,:)
  call gather_all_cr(tmp,nglob,store_val_ux_all,nglob,NPROCTOT)
  ! y-component
  tmp(:) = glob_data(2,:)
  call gather_all_cr(tmp,nglob,store_val_uy_all,nglob,NPROCTOT)
  ! z-component
  tmp(:) = glob_data(3,:)
  call gather_all_cr(tmp,nglob,store_val_uz_all,nglob,NPROCTOT)

  call gather_all_i(ibool,NGLLX*NGLLY*NGLLZ*nspec,ibool_all,NGLLX*NGLLY*NGLLZ*nspec,NPROCTOT)
  call gather_all_i(idoubling,nspec,idoubling_all,nspec,NPROCTOT)


  if (myrank == 0) then

    ! write source and receiver VTK files for Paraview
    open(IOUT_VTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown',iostat=ier)
    if (ier /= 0) stop 'Error opening VTK file'

    write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
    write(IOUT_VTK,'(a)') 'material model VTK file'
    write(IOUT_VTK,'(a)') 'ASCII'
    write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
    write(IOUT_VTK, '(a,i12,a)') 'POINTS ', nglob*NPROCTOT, ' float'
    do iproc = 0, NPROCTOT-1
      do i = 1,nglob

        !x,y,z store have been converted to r theta phi already, need to revert back for xyz output
        rval = store_val_x_all(i,iproc)
        thetaval = store_val_y_all(i,iproc)
        phival = store_val_z_all(i,iproc)
        call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)

        write(IOUT_VTK,'(3e18.6)') xval,yval,zval
      enddo
    enddo
    write(IOUT_VTK,*) ""

    ! defines cell on coarse corner points
    ! note: indices for VTK start at 0
    write(IOUT_VTK,'(a,i12,i12)') "CELLS ",nspec*NPROCTOT,nspec*NPROCTOT*9
    do iproc = 0, NPROCTOT-1
      do ispec = 1,nspec

        ! note: central cube elements are only shared and used in CHUNK_AB and CHUNK_AB_ANTIPODE
        !          all other chunks ignore those elements

        ! specific to inner core elements
        ! exclude fictitious elements in central cube
        if (idoubling_all(ispec,iproc) /= IFLAG_IN_FICTITIOUS_CUBE) then
          ! valid cell
          ! cell corner ids
          write(IOUT_VTK,'(9i12)') 8,ibool_all(1,1,1,ispec,iproc)-1+iproc*nglob, &
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
          write(IOUT_VTK,'(9i12)') 8,ibool_all(1,1,1,1,iproc)-1, &
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
    write(IOUT_VTK,*) ""

    ! type: hexahedrons
    write(IOUT_VTK,'(a,i12)') "CELL_TYPES ",nspec*NPROCTOT
    write(IOUT_VTK,'(6i12)') (12,ispec = 1,nspec*NPROCTOT)
    write(IOUT_VTK,*) ""

    ! x components
    write(IOUT_VTK,'(a,i12)') "POINT_DATA ",nglob*NPROCTOT
    write(IOUT_VTK,'(a)') "SCALARS x_comp float"
    write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
    do iproc = 0, NPROCTOT-1
      do i = 1,nglob
        write(IOUT_VTK,*) store_val_ux_all(i,iproc)
      enddo
    enddo
    ! y components
    write(IOUT_VTK,'(a)') "SCALARS y_comp float"
    write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
    do iproc = 0, NPROCTOT-1
      do i = 1,nglob
        write(IOUT_VTK,*) store_val_uy_all(i,iproc)
      enddo
    enddo
    ! z components
    write(IOUT_VTK,'(a)') "SCALARS z_comp float"
    write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
    do iproc = 0, NPROCTOT-1
      do i = 1,nglob
        write(IOUT_VTK,*) store_val_uz_all(i,iproc)
      enddo
    enddo
    ! norm
    write(IOUT_VTK,'(a)') "SCALARS norm float"
    write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
    do iproc = 0, NPROCTOT-1
      do i = 1,nglob
        write(IOUT_VTK,*) sqrt( store_val_ux_all(i,iproc)**2 &
                          + store_val_uy_all(i,iproc)**2 &
                          + store_val_uz_all(i,iproc)**2 )
      enddo
    enddo
    write(IOUT_VTK,*) ""

    close(IOUT_VTK)

  endif

  deallocate(store_val_x_all,store_val_y_all,store_val_z_all, &
            store_val_ux_all,store_val_uy_all,store_val_uz_all, &
            ibool_all)

  end subroutine write_VTK_data_cr_all

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_VTK_data_elem_cr(nspec,nglob, &
                                    xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                                    gll_data,prname_file)

! external mesh routine for saving vtk files for custom_real values on all gll points

  use constants,only: CUSTOM_REAL,MAX_STRING_LEN,NGLLX,NGLLY,NGLLZ,IOUT_VTK

  implicit none

  integer,intent(in) :: nspec,nglob

  ! global coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob),intent(in) :: xstore_dummy,ystore_dummy,zstore_dummy

  ! gll data values array
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: gll_data

  ! file name
  character(len=MAX_STRING_LEN),intent(in) :: prname_file

  ! local parameters
  integer :: ispec,i,j,k,ier,iglob

  !--------------------------------------------------------------
  ! USER PARAMETERS

  ! flag enabling output on corners only
  logical,parameter :: USE_CORNERS = .false.

  !--------------------------------------------------------------

  !debug
  !print*, '  vtk file: '
  !print*, '    ',prname_file(1:len_trim(prname_file))//'.vtk'

  open(IOUT_VTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown',iostat=ier)
  if (ier /= 0) stop 'Error opening VTK file'

  write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOUT_VTK,'(a)') 'material model VTK file'
  write(IOUT_VTK,'(a)') 'ASCII'
  write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'

  ! writes out all points for each element, not just global ones
  if (USE_CORNERS) then
    write(IOUT_VTK, '(a,i12,a)') 'POINTS ', nspec*8, ' float'
    do ispec=1,nspec
      i = ibool(1,1,1,ispec)
      write(IOUT_VTK,'(3e18.6)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)

      i = ibool(NGLLX,1,1,ispec)
      write(IOUT_VTK,'(3e18.6)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)

      i = ibool(NGLLX,NGLLY,1,ispec)
      write(IOUT_VTK,'(3e18.6)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)

      i = ibool(1,NGLLY,1,ispec)
      write(IOUT_VTK,'(3e18.6)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)

      i = ibool(1,1,NGLLZ,ispec)
      write(IOUT_VTK,'(3e18.6)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)

      i = ibool(NGLLX,1,NGLLZ,ispec)
      write(IOUT_VTK,'(3e18.6)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)

      i = ibool(NGLLX,NGLLY,NGLLZ,ispec)
      write(IOUT_VTK,'(3e18.6)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)

      i = ibool(1,NGLLY,NGLLZ,ispec)
      write(IOUT_VTK,'(3e18.6)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)
    enddo
  else
    write(IOUT_VTK, '(a,i16,a)') 'POINTS ', NGLLX*NGLLY*NGLLZ*nspec, ' float'
    do ispec=1,nspec
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            iglob = ibool(i,j,k,ispec)
            write(IOUT_VTK,'(3e18.6)') xstore_dummy(iglob),ystore_dummy(iglob),zstore_dummy(iglob)
          enddo
        enddo
      enddo
    enddo
  endif
  write(IOUT_VTK,*) ""

  ! note: indices for vtk start at 0
  if (USE_CORNERS) then
    write(IOUT_VTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
    do ispec=1,nspec
      write(IOUT_VTK,'(9i12)') 8, &
            (ispec-1)*8,(ispec-1)*8+1,(ispec-1)*8+2,(ispec-1)*8+3,&
            (ispec-1)*8+4,(ispec-1)*8+5,(ispec-1)*8+6,(ispec-1)*8+7
    enddo
  else
    write(IOUT_VTK,'(a,i16,i16)') "CELLS ",(NGLLX-1)*(NGLLY-1)*(NGLLZ-1)*nspec,(NGLLX-1)*(NGLLY-1)*(NGLLZ-1)*nspec*9
    do ispec=1,nspec
      do k = 1,NGLLZ-1
        do j = 1,NGLLY-1
          do i = 1,NGLLX-1
            write(IOUT_VTK,'(9i16)') 8, &
              (ispec-1)*NGLLZ*NGLLY*NGLLX + (k-1)*NGLLY*NGLLX + (j-1)*NGLLX + i   - 1, &
              (ispec-1)*NGLLZ*NGLLY*NGLLX + (k-1)*NGLLY*NGLLX + (j-1)*NGLLX + i+1 - 1, &
              (ispec-1)*NGLLZ*NGLLY*NGLLX + (k-1)*NGLLY*NGLLX + (j  )*NGLLX + i+1 - 1, &
              (ispec-1)*NGLLZ*NGLLY*NGLLX + (k-1)*NGLLY*NGLLX + (j  )*NGLLX + i   - 1, &
              (ispec-1)*NGLLZ*NGLLY*NGLLX + (k  )*NGLLY*NGLLX + (j-1)*NGLLX + i   - 1, &
              (ispec-1)*NGLLZ*NGLLY*NGLLX + (k  )*NGLLY*NGLLX + (j-1)*NGLLX + i+1 - 1, &
              (ispec-1)*NGLLZ*NGLLY*NGLLX + (k  )*NGLLY*NGLLX + (j  )*NGLLX + i+1 - 1, &
              (ispec-1)*NGLLZ*NGLLY*NGLLX + (k  )*NGLLY*NGLLX + (j  )*NGLLX + i   - 1
          enddo
        enddo
      enddo
    enddo
  endif
  write(IOUT_VTK,*) ""

  ! type: hexahedrons
  if (USE_CORNERS) then
    write(IOUT_VTK,'(a,i12)') "CELL_TYPES ",nspec
    write(IOUT_VTK,'(6i12)') (12,ispec=1,nspec)
  else
    write(IOUT_VTK,'(a,i16)') "CELL_TYPES ",(NGLLX-1)*(NGLLY-1)*(NGLLZ-1)*nspec
    write(IOUT_VTK,'(6i16)') (12,ispec=1,(NGLLX-1)*(NGLLY-1)*(NGLLZ-1)*nspec)
  endif
  write(IOUT_VTK,*) ""

  ! writes out gll-data (velocity) for each element point
  if (USE_CORNERS) then
    write(IOUT_VTK,'(a,i12)') "POINT_DATA ",nspec*8
  else
    write(IOUT_VTK,'(a,i16)') "POINT_DATA ",NGLLX*NGLLY*NGLLZ*nspec
  endif
  write(IOUT_VTK,'(a)') "SCALARS gll_data float"
  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
  if (USE_CORNERS) then
    do ispec = 1,nspec
      !i = ibool(1,1,1,ispec)
      write(IOUT_VTK,*) gll_data(1,1,1,ispec)

      !i = ibool(NGLLX,1,1,ispec)
      write(IOUT_VTK,*) gll_data(NGLLX,1,1,ispec)

      !i = ibool(NGLLX,NGLLY,1,ispec)
      write(IOUT_VTK,*) gll_data(NGLLX,NGLLY,1,ispec)

      !i = ibool(1,NGLLY,1,ispec)
      write(IOUT_VTK,*) gll_data(1,NGLLY,1,ispec)

      !i = ibool(1,1,NGLLZ,ispec)
      write(IOUT_VTK,*) gll_data(1,1,NGLLZ,ispec)

      !i = ibool(NGLLX,1,NGLLZ,ispec)
      write(IOUT_VTK,*) gll_data(NGLLX,1,NGLLZ,ispec)

      !i = ibool(NGLLX,NGLLY,NGLLZ,ispec)
      write(IOUT_VTK,*) gll_data(NGLLX,NGLLY,NGLLZ,ispec)

      !i = ibool(1,NGLLY,NGLLZ,ispec)
      write(IOUT_VTK,*) gll_data(1,NGLLY,NGLLZ,ispec)
    enddo
  else
    do ispec = 1,nspec
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            write(IOUT_VTK,*) gll_data(i,j,k,ispec)
          enddo
        enddo
      enddo
    enddo
  endif
  write(IOUT_VTK,*) ""

  close(IOUT_VTK)

  end subroutine write_VTK_data_elem_cr

