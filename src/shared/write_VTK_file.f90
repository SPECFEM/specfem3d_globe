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

! The VTK writers here are pure Fortran routines, without the need to install and compile the actual
! VTK library files from https://www.vtk.org/download/
!
! VTK encourages to use their library API function calls instead based on an updated VTK installation.
! Still, the routines here are implemented as simple as possible and work well together with Paraview, Visit,
! and other visualization software.
!
! Formats are described here:
! https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
! https://www.vtk.org/Wiki/VTK_XML_Formats#Uncompressed_Data
!
! Two file formats are implemented:
!  .vtk   - the "legacy" VTK format for unstructured grids
!  .vtu   - the "new" VTK XML format for unstructured grids
!
! Most of the routines produce ASCII files, which become much bigger (factor ~4x) than binary file formats.
! For movie files, which will often be generated at multiple time steps, we will use the XML binary format, but
! still produce a single file per time step. This could be improved in future...
!
!
! note: this might be a possible bug in gfortran with -mcmodel=medium on cray,
!       but the write statement
!         write(IOUT_VTK,*) ''
!       produces errors, relocation truncated to fit: R_X86_64_32 against `.lrodata'
!       this can be fixed by using
!         write(IOUT_VTK,*)
!       without quotes to add a newline.


  subroutine write_VTK_data_points(nglob, &
                                   xstore_dummy,ystore_dummy,zstore_dummy, &
                                   points_globalindices,num_points_globalindices, &
                                   prname_file)

! external mesh routine for saving VTK files for points locations

  use constants, only: CUSTOM_REAL,MAX_STRING_LEN,IOUT_VTK

  implicit none

  integer,intent(in) :: nglob

  ! global coordinates
  real(kind=CUSTOM_REAL), dimension(nglob),intent(in) :: xstore_dummy,ystore_dummy,zstore_dummy

  ! GLL data values array
  integer,intent(in) :: num_points_globalindices
  integer, dimension(num_points_globalindices),intent(in) :: points_globalindices

  ! file name
  character(len=MAX_STRING_LEN),intent(in) :: prname_file

  ! local parameters
  integer :: i,iglob,ier

  ! write source and receiver VTK files for Paraview
  !debug
  !write(IMAIN,*) '  VTK file: '
  !write(IMAIN,*) '    ',prname_file(1:len_trim(prname_file))//'.vtk'

  open(IOUT_VTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown',action='write',iostat=ier)
  if (ier /= 0) stop 'Error opening VTK file'

  write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOUT_VTK,'(a)') 'material model VTK file'
  write(IOUT_VTK,'(a)') 'ASCII'
  write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOUT_VTK, '(a,i12,a)') 'POINTS ', num_points_globalindices, ' float'
  do i = 1,num_points_globalindices
    iglob = points_globalindices(i)
    if (iglob <= 0 .or. iglob > nglob) then
      print *,'Error: '//prname_file(1:len_trim(prname_file))//'.vtk'
      print *,'Error global index: ',iglob,i
      stop 'Error VTK points file'
    endif

    write(IOUT_VTK,'(3e18.6)') real(xstore_dummy(iglob),kind=4),real(ystore_dummy(iglob),kind=4),real(zstore_dummy(iglob),kind=4)
  enddo
  write(IOUT_VTK,*)

  close(IOUT_VTK)

  end subroutine write_VTK_data_points

!
!-------------------------------------------------------------------------------------------------
!
!
! used routine, may be used for debugging...
!
!  subroutine write_VTK_glob_points(nglob, &
!                                  xstore_dummy,ystore_dummy,zstore_dummy, &
!                                  glob_values, &
!                                  prname_file)
!
!! external mesh routine for saving VTK files for points locations
!
!  use constants, only: CUSTOM_REAL,MAX_STRING_LEN,IOUT_VTK
!
!  implicit none
!
!  integer,intent(in) :: nglob
!
!  ! global coordinates
!  real(kind=CUSTOM_REAL), dimension(nglob),intent(in) :: xstore_dummy,ystore_dummy,zstore_dummy
!
!  ! GLL data values array
!  real(kind=CUSTOM_REAL), dimension(nglob),intent(in) :: glob_values
!
!  ! file name
!  character(len=MAX_STRING_LEN),intent(in) :: prname_file
!
!  ! local parameters
!  integer :: iglob,ier
!
!  ! write source and receiver VTK files for Paraview
!  !debug
!  !write(IMAIN,*) '  VTK file: '
!  !write(IMAIN,*) '    ',prname_file(1:len_trim(prname_file))//'.vtk'
!
!  open(IOUT_VTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown',action='write',iostat=ier)
!  if (ier /= 0) stop 'Error opening VTK file'
!
!  write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
!  write(IOUT_VTK,'(a)') 'material model VTK file'
!  write(IOUT_VTK,'(a)') 'ASCII'
!  write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
!  write(IOUT_VTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
!  do iglob = 1,nglob
!    write(IOUT_VTK,*) xstore_dummy(iglob),ystore_dummy(iglob),zstore_dummy(iglob)
!  enddo
!  write(IOUT_VTK,*)
!
!  ! writes out GLL data (velocity) for each element point
!  write(IOUT_VTK,'(a,i12)') "POINT_DATA ",nglob
!  write(IOUT_VTK,'(a)') "SCALARS glob_data float"
!  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
!  do iglob = 1,nglob
!    write(IOUT_VTK,*) glob_values(iglob)
!  enddo
!  write(IOUT_VTK,*)
!
!  close(IOUT_VTK)
!
!  end subroutine write_VTK_glob_points
!
!
!
!-------------------------------------------------------------------------------------------------
!


  subroutine write_VTK_data_elem_l(nspec,nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                        elem_flag,prname_file)

! routine for saving VTK file holding logical flag on each spectral element

  use constants, only: CUSTOM_REAL,MAX_STRING_LEN,NGLLX,NGLLY,NGLLZ,IOUT_VTK

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

  open(IOUT_VTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown',action='write',iostat=ier)
  if (ier /= 0) stop 'Error opening VTK file'

  write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOUT_VTK,'(a)') 'material model VTK file'
  write(IOUT_VTK,'(a)') 'ASCII'
  write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOUT_VTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do i = 1,nglob
    write(IOUT_VTK,'(3e18.6)') real(xstore_dummy(i),kind=4),real(ystore_dummy(i),kind=4),real(zstore_dummy(i),kind=4)
  enddo
  write(IOUT_VTK,*)

  ! note: indices for VTK start at 0
  write(IOUT_VTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
  do ispec = 1,nspec
    write(IOUT_VTK,'(9i12)') 8, &
          ibool(1,1,1,ispec)-1,ibool(NGLLX,1,1,ispec)-1,ibool(NGLLX,NGLLY,1,ispec)-1,ibool(1,NGLLY,1,ispec)-1, &
          ibool(1,1,NGLLZ,ispec)-1,ibool(NGLLX,1,NGLLZ,ispec)-1,ibool(NGLLX,NGLLY,NGLLZ,ispec)-1,ibool(1,NGLLY,NGLLZ,ispec)-1
  enddo
  write(IOUT_VTK,*)

  ! type: hexahedrons
  write(IOUT_VTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOUT_VTK,'(6i12)') (12,ispec = 1,nspec)
  write(IOUT_VTK,*)

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
  write(IOUT_VTK,*)
  close(IOUT_VTK)

  end subroutine write_VTK_data_elem_l


!
!-------------------------------------------------------------------------------------------------
!


  subroutine write_VTK_data_elem_i(nspec,nglob, &
                                   xstore_dummy,ystore_dummy,zstore_dummy, &
                                   ibool,elem_flag,prname_file)


! routine for saving VTK file holding integer value on each spectral element

  use constants, only: CUSTOM_REAL,MAX_STRING_LEN,NGLLX,NGLLY,NGLLZ,IOUT_VTK

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

  open(IOUT_VTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown',action='write',iostat=ier)
  if (ier /= 0) stop 'Error opening VTK file'

  write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOUT_VTK,'(a)') 'material model VTK file'
  write(IOUT_VTK,'(a)') 'ASCII'
  write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOUT_VTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do i = 1,nglob
    write(IOUT_VTK,'(3e18.6)') real(xstore_dummy(i),kind=4),real(ystore_dummy(i),kind=4),real(zstore_dummy(i),kind=4)
  enddo
  write(IOUT_VTK,*)

  ! note: indices for VTK start at 0
  write(IOUT_VTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
  do ispec = 1,nspec
    write(IOUT_VTK,'(9i12)') 8, &
          ibool(1,1,1,ispec)-1,ibool(NGLLX,1,1,ispec)-1,ibool(NGLLX,NGLLY,1,ispec)-1,ibool(1,NGLLY,1,ispec)-1, &
          ibool(1,1,NGLLZ,ispec)-1,ibool(NGLLX,1,NGLLZ,ispec)-1,ibool(NGLLX,NGLLY,NGLLZ,ispec)-1,ibool(1,NGLLY,NGLLZ,ispec)-1
  enddo
  write(IOUT_VTK,*)

  ! type: hexahedrons
  write(IOUT_VTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOUT_VTK,'(6i12)') (12,ispec = 1,nspec)
  write(IOUT_VTK,*)

  write(IOUT_VTK,'(a,i12)') "CELL_DATA ",nspec
  write(IOUT_VTK,'(a)') "SCALARS elem_val integer"
  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
  do ispec = 1,nspec
    write(IOUT_VTK,*) elem_flag(ispec)
  enddo
  write(IOUT_VTK,*)
  close(IOUT_VTK)

  end subroutine write_VTK_data_elem_i

!
!-------------------------------------------------------------------------------------------------
!

! external mesh routine for saving VTK files for CUSTOM_REAL values on global points

  subroutine write_VTK_data_cr(idoubling,nspec,nglob, &
                               rstore_dummy, &
                               ibool,glob_data,prname_file)

! outputs single file for each process

  use constants, only: CUSTOM_REAL,MAX_STRING_LEN,NDIM,NGLLX,NGLLY,NGLLZ,IOUT_VTK,IFLAG_IN_FICTITIOUS_CUBE

  implicit none

  integer,intent(in) :: nspec,nglob

  integer, dimension(nspec),intent(in) :: idoubling

  ! global coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(in) :: rstore_dummy

  ! global data values array
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(in) :: glob_data

  ! file name
  character(len=MAX_STRING_LEN),intent(in) :: prname_file

  ! local parameters
  integer :: ispec,i,ier
  real(kind=CUSTOM_REAL) :: rval,thetaval,phival,xval,yval,zval

  ! write source and receiver VTK files for Paraview
  open(IOUT_VTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown',action='write',iostat=ier)
  if (ier /= 0) stop 'Error opening VTK file'

  write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOUT_VTK,'(a)') 'material model VTK file'
  write(IOUT_VTK,'(a)') 'ASCII'
  write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOUT_VTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do i = 1,nglob

    !x,y,z store have been converted to r theta phi already, need to revert back for xyz output
    rval = rstore_dummy(1,i)
    thetaval = rstore_dummy(2,i)
    phival = rstore_dummy(3,i)

    call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)

    write(IOUT_VTK,'(3e18.6)') real(xval,kind=4),real(yval,kind=4),real(zval,kind=4)
  enddo
  write(IOUT_VTK,*)

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
  write(IOUT_VTK,*)

  ! type: hexahedrons
  write(IOUT_VTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOUT_VTK,'(6i12)') (12,ispec = 1,nspec)
  write(IOUT_VTK,*)

  ! x components
  write(IOUT_VTK,'(a,i12)') "POINT_DATA ",nglob
  write(IOUT_VTK,'(a)') "SCALARS x_comp float"
  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
  do i = 1,nglob
      write(IOUT_VTK,*) real(glob_data(1,i),kind=4)
  enddo
  ! y components
  write(IOUT_VTK,'(a)') "SCALARS y_comp float"
  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
  do i = 1,nglob
      write(IOUT_VTK,*) real(glob_data(2,i),kind=4)
  enddo
  ! z components
  write(IOUT_VTK,'(a)') "SCALARS z_comp float"
  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
  do i = 1,nglob
      write(IOUT_VTK,*) real(glob_data(3,i),kind=4)
  enddo
  ! norm
  write(IOUT_VTK,'(a)') "SCALARS norm float"
  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
  do i = 1,nglob
      write(IOUT_VTK,*) sqrt( real(glob_data(1,i)*glob_data(1,i) &
                                 + glob_data(2,i)*glob_data(2,i) &
                                 + glob_data(3,i)*glob_data(3,i),kind=4))
  enddo
  write(IOUT_VTK,*)

  close(IOUT_VTK)


  end subroutine write_VTK_data_cr

!
!-------------------------------------------------------------------------------------------------
!
!
! unused routine, may be used for debugging...
!
!! external mesh routine for saving VTK files for CUSTOM_REAL values on global points
!
!  subroutine write_VTK_data_cr_all(myrank,NPROCTOT,idoubling, &
!                              nspec,nglob, &
!                              xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
!                              glob_data,prname_file)
!
!! outputs single file for all processes
!
!  use constants, only: CUSTOM_REAL,MAX_STRING_LEN,NDIM,NGLLX,NGLLY,NGLLZ,IOUT_VTK,IFLAG_IN_FICTITIOUS_CUBE
!
!  implicit none
!
!  integer,intent(in) :: myrank,NPROCTOT
!
!  integer,intent(in) ::nspec,nglob
!
!  integer, dimension(nspec),intent(in) :: idoubling
!
!  ! global coordinates
!  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool
!  real(kind=CUSTOM_REAL), dimension(nglob),intent(in) :: xstore_dummy,ystore_dummy,zstore_dummy
!
!  ! global data values array
!  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(in) :: glob_data
!
!  ! file name
!  character(len=MAX_STRING_LEN),intent(in) :: prname_file
!
!  ! local parameters
!  integer :: ispec,i,iproc,ier
!  real(kind=CUSTOM_REAL) :: rval,thetaval,phival,xval,yval,zval
!
!  real(kind=CUSTOM_REAL), dimension(:,:),allocatable :: &
!      store_val_x_all,store_val_y_all,store_val_z_all, &
!      store_val_ux_all,store_val_uy_all,store_val_uz_all
!  integer, dimension(:,:,:,:,:),allocatable :: ibool_all
!  integer, dimension(:,:),allocatable :: idoubling_all
!  real(kind=CUSTOM_REAL), dimension(nglob) :: tmp
!
!  ! main collect arrays
!  if (myrank == 0) then
!    allocate(store_val_x_all(nglob,0:NPROCTOT-1), &
!            store_val_y_all(nglob,0:NPROCTOT-1), &
!            store_val_z_all(nglob,0:NPROCTOT-1), &
!            store_val_ux_all(nglob,0:NPROCTOT-1), &
!            store_val_uy_all(nglob,0:NPROCTOT-1), &
!            store_val_uz_all(nglob,0:NPROCTOT-1), &
!            idoubling_all(nspec,0:NPROCTOT-1), &
!            ibool_all(NGLLX,NGLLY,NGLLZ,nspec,0:NPROCTOT-1),stat=ier)
!    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating stores')
!  else
!    ! dummy arrays
!    allocate(store_val_x_all(1,1), &
!            store_val_y_all(1,1), &
!            store_val_z_all(1,1), &
!            store_val_ux_all(1,1), &
!            store_val_uy_all(1,1), &
!            store_val_uz_all(1,1), &
!            idoubling_all(1,1), &
!            ibool_all(1,1,1,1,1),stat=ier)
!    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating dummy stores')
!  endif
!
!  ! gather info on main proc
!  call gather_all_cr(xstore_dummy,nglob,store_val_x_all,nglob,NPROCTOT)
!  call gather_all_cr(ystore_dummy,nglob,store_val_y_all,nglob,NPROCTOT)
!  call gather_all_cr(zstore_dummy,nglob,store_val_z_all,nglob,NPROCTOT)
!
!  ! attention: these calls produce copies of the glob_data array
!  ! x-component
!  tmp(:) = glob_data(1,:)
!  call gather_all_cr(tmp,nglob,store_val_ux_all,nglob,NPROCTOT)
!  ! y-component
!  tmp(:) = glob_data(2,:)
!  call gather_all_cr(tmp,nglob,store_val_uy_all,nglob,NPROCTOT)
!  ! z-component
!  tmp(:) = glob_data(3,:)
!  call gather_all_cr(tmp,nglob,store_val_uz_all,nglob,NPROCTOT)
!
!  call gather_all_i(ibool,NGLLX*NGLLY*NGLLZ*nspec,ibool_all,NGLLX*NGLLY*NGLLZ*nspec,NPROCTOT)
!  call gather_all_i(idoubling,nspec,idoubling_all,nspec,NPROCTOT)
!
!
!  if (myrank == 0) then
!
!    ! write source and receiver VTK files for Paraview
!    open(IOUT_VTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown',action='write',iostat=ier)
!    if (ier /= 0) stop 'Error opening VTK file'
!
!    write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
!    write(IOUT_VTK,'(a)') 'material model VTK file'
!    write(IOUT_VTK,'(a)') 'ASCII'
!    write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
!    write(IOUT_VTK, '(a,i12,a)') 'POINTS ', nglob*NPROCTOT, ' float'
!    do iproc = 0, NPROCTOT-1
!      do i = 1,nglob
!
!        !x,y,z store have been converted to r theta phi already, need to revert back for xyz output
!        rval = store_val_x_all(i,iproc)
!        thetaval = store_val_y_all(i,iproc)
!        phival = store_val_z_all(i,iproc)
!        call rthetaphi_2_xyz(xval,yval,zval,rval,thetaval,phival)
!
!        write(IOUT_VTK,'(3e18.6)') xval,yval,zval
!      enddo
!    enddo
!    write(IOUT_VTK,*)
!
!    ! defines cell on coarse corner points
!    ! note: indices for VTK start at 0
!    write(IOUT_VTK,'(a,i12,i12)') "CELLS ",nspec*NPROCTOT,nspec*NPROCTOT*9
!    do iproc = 0, NPROCTOT-1
!      do ispec = 1,nspec
!
!        ! note: central cube elements are only shared and used in CHUNK_AB and CHUNK_AB_ANTIPODE
!        !          all other chunks ignore those elements
!
!        ! specific to inner core elements
!        ! exclude fictitious elements in central cube
!        if (idoubling_all(ispec,iproc) /= IFLAG_IN_FICTITIOUS_CUBE) then
!          ! valid cell
!          ! cell corner ids
!          write(IOUT_VTK,'(9i12)') 8,ibool_all(1,1,1,ispec,iproc)-1+iproc*nglob, &
!                            ibool_all(NGLLX,1,1,ispec,iproc)-1+iproc*nglob, &
!                            ibool_all(NGLLX,NGLLY,1,ispec,iproc)-1+iproc*nglob, &
!                            ibool_all(1,NGLLY,1,ispec,iproc)-1+iproc*nglob, &
!                            ibool_all(1,1,NGLLZ,ispec,iproc)-1+iproc*nglob, &
!                            ibool_all(NGLLX,1,NGLLZ,ispec,iproc)-1+iproc*nglob, &
!                            ibool_all(NGLLX,NGLLY,NGLLZ,ispec,iproc)-1+iproc*nglob, &
!                            ibool_all(1,NGLLY,NGLLZ,ispec,iproc)-1+iproc*nglob
!        else
!          ! fictitious elements in central cube
!          ! maps cell onto a randomly chosen point
!          write(IOUT_VTK,'(9i12)') 8,ibool_all(1,1,1,1,iproc)-1, &
!                            ibool_all(1,1,1,1,iproc)-1, &
!                            ibool_all(1,1,1,1,iproc)-1, &
!                            ibool_all(1,1,1,1,iproc)-1, &
!                            ibool_all(1,1,1,1,iproc)-1, &
!                            ibool_all(1,1,1,1,iproc)-1, &
!                            ibool_all(1,1,1,1,iproc)-1, &
!                            ibool_all(1,1,1,1,iproc)-1
!        endif
!
!      enddo
!    enddo
!    write(IOUT_VTK,*)
!
!    ! type: hexahedrons
!    write(IOUT_VTK,'(a,i12)') "CELL_TYPES ",nspec*NPROCTOT
!    write(IOUT_VTK,'(6i12)') (12,ispec = 1,nspec*NPROCTOT)
!    write(IOUT_VTK,*)
!
!    ! x components
!    write(IOUT_VTK,'(a,i12)') "POINT_DATA ",nglob*NPROCTOT
!    write(IOUT_VTK,'(a)') "SCALARS x_comp float"
!    write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
!    do iproc = 0, NPROCTOT-1
!      do i = 1,nglob
!        write(IOUT_VTK,*) store_val_ux_all(i,iproc)
!      enddo
!    enddo
!    ! y components
!    write(IOUT_VTK,'(a)') "SCALARS y_comp float"
!    write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
!    do iproc = 0, NPROCTOT-1
!      do i = 1,nglob
!        write(IOUT_VTK,*) store_val_uy_all(i,iproc)
!      enddo
!    enddo
!    ! z components
!    write(IOUT_VTK,'(a)') "SCALARS z_comp float"
!    write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
!    do iproc = 0, NPROCTOT-1
!      do i = 1,nglob
!        write(IOUT_VTK,*) store_val_uz_all(i,iproc)
!      enddo
!    enddo
!    ! norm
!    write(IOUT_VTK,'(a)') "SCALARS norm float"
!    write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
!    do iproc = 0, NPROCTOT-1
!      do i = 1,nglob
!        write(IOUT_VTK,*) sqrt( store_val_ux_all(i,iproc)**2 &
!                          + store_val_uy_all(i,iproc)**2 &
!                          + store_val_uz_all(i,iproc)**2 )
!      enddo
!    enddo
!    write(IOUT_VTK,*)
!
!    close(IOUT_VTK)
!
!  endif
!
!  deallocate(store_val_x_all,store_val_y_all,store_val_z_all, &
!            store_val_ux_all,store_val_uy_all,store_val_uz_all, &
!            ibool_all)
!
!  end subroutine write_VTK_data_cr_all
!
!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_VTK_data_gll_cr(nspec,nglob, &
                                    xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                                    gll_data,prname_file)

! external mesh routine for saving vtk files for CUSTOM_REAL values on all GLL points

  use constants, only: CUSTOM_REAL,MAX_STRING_LEN,NGLLX,NGLLY,NGLLZ,IOUT_VTK

  implicit none

  integer,intent(in) :: nspec,nglob

  ! global coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob),intent(in) :: xstore_dummy,ystore_dummy,zstore_dummy

  ! GLL data values array
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
  !print *, '  vtk file: '
  !print *, '    ',prname_file(1:len_trim(prname_file))//'.vtk'

  open(IOUT_VTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown',action='write',iostat=ier)
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
      write(IOUT_VTK,'(3e18.6)') real(xstore_dummy(i),kind=4),real(ystore_dummy(i),kind=4),real(zstore_dummy(i),kind=4)

      i = ibool(NGLLX,1,1,ispec)
      write(IOUT_VTK,'(3e18.6)') real(xstore_dummy(i),kind=4),real(ystore_dummy(i),kind=4),real(zstore_dummy(i),kind=4)

      i = ibool(NGLLX,NGLLY,1,ispec)
      write(IOUT_VTK,'(3e18.6)') real(xstore_dummy(i),kind=4),real(ystore_dummy(i),kind=4),real(zstore_dummy(i),kind=4)

      i = ibool(1,NGLLY,1,ispec)
      write(IOUT_VTK,'(3e18.6)') real(xstore_dummy(i),kind=4),real(ystore_dummy(i),kind=4),real(zstore_dummy(i),kind=4)

      i = ibool(1,1,NGLLZ,ispec)
      write(IOUT_VTK,'(3e18.6)') real(xstore_dummy(i),kind=4),real(ystore_dummy(i),kind=4),real(zstore_dummy(i),kind=4)

      i = ibool(NGLLX,1,NGLLZ,ispec)
      write(IOUT_VTK,'(3e18.6)') real(xstore_dummy(i),kind=4),real(ystore_dummy(i),kind=4),real(zstore_dummy(i),kind=4)

      i = ibool(NGLLX,NGLLY,NGLLZ,ispec)
      write(IOUT_VTK,'(3e18.6)') real(xstore_dummy(i),kind=4),real(ystore_dummy(i),kind=4),real(zstore_dummy(i),kind=4)

      i = ibool(1,NGLLY,NGLLZ,ispec)
      write(IOUT_VTK,'(3e18.6)') real(xstore_dummy(i),kind=4),real(ystore_dummy(i),kind=4),real(zstore_dummy(i),kind=4)
    enddo
  else
    write(IOUT_VTK, '(a,i16,a)') 'POINTS ', NGLLX*NGLLY*NGLLZ*nspec, ' float'
    do ispec=1,nspec
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            iglob = ibool(i,j,k,ispec)
            write(IOUT_VTK,'(3e18.6)') real(xstore_dummy(iglob),kind=4), &
                                       real(ystore_dummy(iglob),kind=4), &
                                       real(zstore_dummy(iglob),kind=4)
          enddo
        enddo
      enddo
    enddo
  endif
  write(IOUT_VTK,*)

  ! note: indices for vtk start at 0
  if (USE_CORNERS) then
    write(IOUT_VTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
    do ispec=1,nspec
      write(IOUT_VTK,'(9i12)') 8, &
            (ispec-1)*8,(ispec-1)*8+1,(ispec-1)*8+2,(ispec-1)*8+3, &
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
  write(IOUT_VTK,*)

  ! type: hexahedrons
  if (USE_CORNERS) then
    write(IOUT_VTK,'(a,i12)') "CELL_TYPES ",nspec
    write(IOUT_VTK,'(6i12)') (12,ispec=1,nspec)
  else
    write(IOUT_VTK,'(a,i16)') "CELL_TYPES ",(NGLLX-1)*(NGLLY-1)*(NGLLZ-1)*nspec
    write(IOUT_VTK,'(6i16)') (12,ispec=1,(NGLLX-1)*(NGLLY-1)*(NGLLZ-1)*nspec)
  endif
  write(IOUT_VTK,*)

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
      write(IOUT_VTK,*) real(gll_data(1,1,1,ispec),kind=4)

      !i = ibool(NGLLX,1,1,ispec)
      write(IOUT_VTK,*) real(gll_data(NGLLX,1,1,ispec),kind=4)

      !i = ibool(NGLLX,NGLLY,1,ispec)
      write(IOUT_VTK,*) real(gll_data(NGLLX,NGLLY,1,ispec),kind=4)

      !i = ibool(1,NGLLY,1,ispec)
      write(IOUT_VTK,*) real(gll_data(1,NGLLY,1,ispec),kind=4)

      !i = ibool(1,1,NGLLZ,ispec)
      write(IOUT_VTK,*) real(gll_data(1,1,NGLLZ,ispec),kind=4)

      !i = ibool(NGLLX,1,NGLLZ,ispec)
      write(IOUT_VTK,*) real(gll_data(NGLLX,1,NGLLZ,ispec),kind=4)

      !i = ibool(NGLLX,NGLLY,NGLLZ,ispec)
      write(IOUT_VTK,*) real(gll_data(NGLLX,NGLLY,NGLLZ,ispec),kind=4)

      !i = ibool(1,NGLLY,NGLLZ,ispec)
      write(IOUT_VTK,*) real(gll_data(1,NGLLY,NGLLZ,ispec),kind=4)
    enddo
  else
    do ispec = 1,nspec
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            write(IOUT_VTK,*) real(gll_data(i,j,k,ispec),kind=4)
          enddo
        enddo
      enddo
    enddo
  endif
  write(IOUT_VTK,*)

  close(IOUT_VTK)

  end subroutine write_VTK_data_gll_cr


!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_VTK_data_elem_cr(nspec,nglob, &
                                    xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                                    elem_data,prname_file)

! saves vtk files for CUSTOM_REAL values on all spectral elements

  use constants, only: CUSTOM_REAL,MAX_STRING_LEN,NGLLX,NGLLY,NGLLZ,IOUT_VTK

  implicit none

  integer,intent(in) :: nspec,nglob

  ! global coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob),intent(in) :: xstore_dummy,ystore_dummy,zstore_dummy

  ! element data values array
  real(kind=CUSTOM_REAL), dimension(nspec),intent(in) :: elem_data

  ! file name
  character(len=MAX_STRING_LEN),intent(in) :: prname_file

  ! local parameters
  integer :: ispec,i,ier

  open(IOUT_VTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown',action='write',iostat=ier)
  if (ier /= 0) stop 'Error opening VTK file'

  write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOUT_VTK,'(a)') 'material model VTK file'
  write(IOUT_VTK,'(a)') 'ASCII'
  write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOUT_VTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do i=1,nglob
    write(IOUT_VTK,'(3e18.6)') real(xstore_dummy(i),kind=4),real(ystore_dummy(i),kind=4),real(zstore_dummy(i),kind=4)
  enddo
  write(IOUT_VTK,*) ''

  ! note: indices for vtk start at 0
  write(IOUT_VTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
  do ispec = 1,nspec
    write(IOUT_VTK,'(9i12)') 8, &
          ibool(1,1,1,ispec)-1,ibool(NGLLX,1,1,ispec)-1,ibool(NGLLX,NGLLY,1,ispec)-1,ibool(1,NGLLY,1,ispec)-1, &
          ibool(1,1,NGLLZ,ispec)-1,ibool(NGLLX,1,NGLLZ,ispec)-1,ibool(NGLLX,NGLLY,NGLLZ,ispec)-1,ibool(1,NGLLY,NGLLZ,ispec)-1
  enddo
  write(IOUT_VTK,*) ''

  ! type: hexahedrons
  write(IOUT_VTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOUT_VTK,'(6i12)') (12,ispec=1,nspec)
  write(IOUT_VTK,*) ''

  write(IOUT_VTK,'(a,i12)') "CELL_DATA ",nspec
  write(IOUT_VTK,'(a)') "SCALARS elem_val float"
  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
  do ispec = 1,nspec
    write(IOUT_VTK,*) real(elem_data(ispec),kind=4)
  enddo
  write(IOUT_VTK,*) ''
  close(IOUT_VTK)

  end subroutine write_VTK_data_elem_cr


!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_VTU_data_elem_cr_binary(nspec,nglob, &
                                           xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                                           elem_data,prname_file)

! saves vtu files in binary format, for CUSTOM_REAL values on all spectral elements

  use constants, only: CUSTOM_REAL,SIZE_REAL,SIZE_INTEGER,MAX_STRING_LEN,NGLLX,NGLLY,NGLLZ,IOUT_VTK

  implicit none

  integer,intent(in) :: nspec,nglob

  ! global coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob),intent(in) :: xstore_dummy,ystore_dummy,zstore_dummy

  ! element data values array
  real(kind=CUSTOM_REAL), dimension(nspec),intent(in) :: elem_data

  ! file name
  character(len=MAX_STRING_LEN),intent(in) :: prname_file

  ! local parameters
  integer :: ispec,i,it,ier

  ! local parameters
  integer :: len_bytes,offset
  character(len=MAX_STRING_LEN) :: var_name
  character(len=16) :: str1,str_endian
  character(len=12) :: str2,str3
  character(len=1),parameter :: LF = achar(10) ! line feed character

  ! endianness - determine endianness at run time:
  ! https://www.ibm.com/developerworks/aix/library/au-endianc/index.html
  !
  ! for the Fortran version:
  ! given the integer value of 1, big endian uses 00 00 00 01 and little endian uses 01 00 00 00 as bit representation.
  ! using transfer(1,'a') takes the bit-representation of integer value 1 (multi-byte, usually 32-bit)
  ! and interprets it as a character type (of 1-byte length).
  ! thus, looking at the first byte, it is either 0 or 1, respectively.
  ! finally, ichar(c) takes a character and returns its position number.
  logical, parameter :: is_big_endian = (ichar(transfer(1,'a')) == 0)

  ! note: VTK by default assumes binary data is big endian for "legacy" VTK format,
  !       we use here the new .vtu file with XML format. in this case, we can specify the endianness by byte_order=".."

  if (is_big_endian) then
    str_endian = 'BigEndian'
  else
    str_endian = 'LittleEndian'
  endif

  ! variable name
  var_name = 'elem_val'

  ! numbers as strings
  write(str1,'(i16)') nglob
  write(str2,'(i12)') nspec

  ! data offset for appended datasets
  offset = 0
  write(str3,'(i12)') offset

  ! opens unstructured grid file as binary file
  ! convert='BIG_ENDIAN' not needed, will be specified in XML format
  open(IOUT_VTK,file=trim(prname_file)//'.vtu',access='stream',form='unformatted', &
         status='unknown',action='write',iostat=ier)
  if (ier /= 0 ) then
    print *, 'Error opening VTU output file: ',trim(prname_file)//'.vtu'
    stop 'Error opening VTU output file'
  endif

  ! XML file header
  ! see: https://www.vtk.org/Wiki/VTK_XML_Formats
  write(IOUT_VTK) '<?xml version="1.0"?>' // LF
  write(IOUT_VTK) '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="'// trim(str_endian) // '">' // LF
  write(IOUT_VTK) '<UnstructuredGrid>' // LF
  write(IOUT_VTK) '<Piece NumberOfPoints="'// str1 // '" NumberOfCells="' // str2 // '">' // LF

  ! points
  write(IOUT_VTK) '<Points>' // LF
  ! binary format: not working properly yet - using appended instead
  !write(IOUT_VTK) '<DataArray type="Float32" Name="Points" NumberOfComponents="3" format="binary" encoding="raw">' // LF
  !! number of bytes to follow
  !! see format description: https://www.vtk.org/Wiki/VTK_XML_Formats#Uncompressed_Data
  !len_bytes = 3 * nglob * SIZE_REAL
  !write(IOUT_VTK) len_bytes
  !do i = 1,nglob
  !  write(IOUT_VTK) real(total_dat_xyz(1,i),kind=4),real(total_dat_xyz(2,i),kind=4),real(total_dat_xyz(3,i),kind=4)
  !enddo
  !write(IOUT_VTK) '</DataArray>' // LF
  !
  ! appended format:
  write(IOUT_VTK) '<DataArray type="Float32" Name="Points" NumberOfComponents="3" format="appended" offset="' &
                   // str3 // '"/>' // LF
  ! updates offset
  ! array length in bytes
  len_bytes = 3 * nglob * SIZE_REAL
  ! new offset: data array size (len_bytes) + 4 bytes for length specifier at the beginning
  offset = offset + len_bytes + 4
  write(str3,'(i12)') offset
  write(IOUT_VTK) '</Points>'//LF

  ! cells
  write(IOUT_VTK) '<Cells>' // LF
  ! connectivity
  write(IOUT_VTK) '<DataArray type="Int32" Name="connectivity" format="appended" offset="' // str3 // '"/>' // LF
  ! updates offset
  len_bytes = 8 * nspec * SIZE_INTEGER
  offset = offset + len_bytes + 4
  write(str3,'(i12)') offset

  ! offsets
  write(IOUT_VTK) '<DataArray type="Int32" Name="offsets" format="appended" offset="' // str3 // '"/>' // LF
  ! updates offset
  len_bytes = nspec * SIZE_INTEGER
  offset = offset + len_bytes + 4
  write(str3,'(i12)') offset

  ! type: hexahedrons
  write(IOUT_VTK) '<DataArray type="Int32" Name="types" format="appended" offset="' // str3 // '"/>' // LF
  ! updates offset
  len_bytes = nspec * SIZE_INTEGER
  offset = offset + len_bytes + 4
  write(str3,'(i12)') offset
  write(IOUT_VTK) '</Cells>' // LF

  ! cell data
  write(IOUT_VTK) '<CellData Scalars="Scalars_">' // LF
  write(IOUT_VTK) '<DataArray type="Float32" Name="' // trim(var_name) // '" format="appended" offset="' &
                  // str3 // '"/>' // LF
  write(IOUT_VTK) '</CellData>' // LF

  ! empty point data values
  write(IOUT_VTK) '<PointData>' // '</PointData>' // LF

  ! finishes XML file
  write(IOUT_VTK) '</Piece>' // LF
  write(IOUT_VTK) '</UnstructuredGrid>' // LF

  ! in case of appended data format, with offsets:
  !write(IOUT_VTK) '<AppendedData encoding="base64">' // LF
  write(IOUT_VTK) '<AppendedData encoding="raw">' // LF
  write(IOUT_VTK) '_'
  ! points
  len_bytes = 3 * nglob * SIZE_REAL
  write(IOUT_VTK) len_bytes
  do i = 1,nglob
    write(IOUT_VTK) real(xstore_dummy(i),kind=4),real(ystore_dummy(i),kind=4),real(zstore_dummy(i),kind=4)
  enddo
  ! cells
  ! connectivity
  ! number of bytes to follow
  len_bytes = 8 * nspec * SIZE_INTEGER
  write(IOUT_VTK) len_bytes
  ! note: indices for VTK start at 0
  do ispec = 1,nspec
    write(IOUT_VTK) ibool(1,1,1,ispec)-1,ibool(NGLLX,1,1,ispec)-1,ibool(NGLLX,NGLLY,1,ispec)-1,ibool(1,NGLLY,1,ispec)-1, &
      ibool(1,1,NGLLZ,ispec)-1,ibool(NGLLX,1,NGLLZ,ispec)-1,ibool(NGLLX,NGLLY,NGLLZ,ispec)-1,ibool(1,NGLLY,NGLLZ,ispec)-1
  enddo
  ! offsets
  ! number of bytes to follow
  len_bytes = nspec * SIZE_INTEGER
  write(IOUT_VTK) len_bytes
  ! offsets (8 points each)
  write(IOUT_VTK) ((it*8),it = 1,nspec)
  ! types
  ! number of bytes to follow
  len_bytes = nspec * SIZE_INTEGER
  write(IOUT_VTK) len_bytes
  ! type for hex elements is 12
  write(IOUT_VTK) (12,it = 1,nspec)

  ! cell data values
  ! number of bytes to follow
  len_bytes = nspec * SIZE_REAL
  write(IOUT_VTK) len_bytes
  ! data values
  do ispec = 1,nspec
    write(IOUT_VTK) real(elem_data(ispec),kind=4)
  enddo

  write(IOUT_VTK) '</AppendedData>' // LF
  write(IOUT_VTK) '</VTKFile>' // LF
  close(IOUT_VTK)

  end subroutine write_VTU_data_elem_cr_binary


!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_VTK_movie_data(ne,np,total_dat_xyz,total_dat_con,total_dat,mesh_file,var_name)

! saves vtk file for CUSTOM_REAL values

  use constants, only: CUSTOM_REAL,MAX_STRING_LEN,IOUT_VTK,SIZE_DOUBLE

  implicit none

  integer,intent(in) :: ne,np

  ! coordinates
  real(kind=CUSTOM_REAL), dimension(3,np),intent(in) :: total_dat_xyz
  ! connections
  integer, dimension(8,ne),intent(in) :: total_dat_con
  ! data values array
  real(kind=CUSTOM_REAL), dimension(np),intent(in) :: total_dat
  ! file name
  character(len=MAX_STRING_LEN),intent(in) :: mesh_file,var_name

  ! local parameters
  integer :: i,it,ier
  real :: val

  ! opens unstructured grid file
  open(IOUT_VTK,file=mesh_file(1:len_trim(mesh_file)),status='unknown',iostat=ier)
  if (ier /= 0 ) then
    print *, 'Error opening VTK output file: ',trim(mesh_file)
    stop 'Error opening VTK output file'
  endif
  write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOUT_VTK,'(a)') 'material model VTK file'
  write(IOUT_VTK,'(a)') 'ASCII'
  write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'

  ! points
  write(IOUT_VTK, '(a,i16,a)') 'POINTS ', np, ' float'
  do i = 1,np
    write(IOUT_VTK,'(3e18.6)') real(total_dat_xyz(1,i),kind=4),real(total_dat_xyz(2,i),kind=4),real(total_dat_xyz(3,i),kind=4)
  enddo
  write(IOUT_VTK,*) ''

  ! cells
  ! note: indices for VTK start at 0
  write(IOUT_VTK,'(a,i12,i12)') "CELLS ",ne,ne*9
  do i = 1,ne
    write(IOUT_VTK,'(9i12)') 8,total_dat_con(1,i),total_dat_con(2,i),total_dat_con(3,i),total_dat_con(4,i), &
                               total_dat_con(5,i),total_dat_con(6,i),total_dat_con(7,i),total_dat_con(8,i)
  enddo
  write(IOUT_VTK,*) ''

  ! type: hexahedrons
  write(IOUT_VTK,'(a,i12)') "CELL_TYPES ",ne
  write(IOUT_VTK,'(6i12)') (12,it = 1,ne)
  write(IOUT_VTK,*) ''

  ! data values
  write(IOUT_VTK,'(a,i12)') "POINT_DATA ",np
  write(IOUT_VTK,'(a)') "SCALARS "//trim(var_name)//" float"
  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
  if (CUSTOM_REAL == SIZE_DOUBLE) then
    ! double precision values
    do i = 1,np
      ! converts to float
      val = real(total_dat(i),kind=4)
      ! stay within boundaries of float values, otherwise paraview will complain
      if (abs(val) < tiny(val)) val = sign(1.0,val) * tiny(val)
      if (abs(val) > huge(val)) val = sign(1.0,val) * huge(val)
      write(IOUT_VTK,*) val
    enddo
  else
    ! single precision
    do i = 1,np
      write(IOUT_VTK,*) total_dat(i)
    enddo
  endif
  write(IOUT_VTK,*) ''
  close(IOUT_VTK)

  end subroutine write_VTK_movie_data

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_VTK_movie_data_elemental(ne,np,total_dat_xyz,total_dat_con,total_dat,mesh_file,var_name)

! saves vtk file for CUSTOM_REAL values, writes out each element as unconnected finite elements
! (no shared global points). this visualizes sharp jumps/discontinuities from one element to another.

  use constants, only: CUSTOM_REAL,MAX_STRING_LEN,IOUT_VTK,SIZE_DOUBLE

  implicit none

  integer,intent(in) :: ne,np

  ! coordinates
  real(kind=CUSTOM_REAL), dimension(3,np),intent(in) :: total_dat_xyz
  ! connections
  integer, dimension(8,ne),intent(in) :: total_dat_con
  ! data values array
  real(kind=CUSTOM_REAL), dimension(np),intent(in) :: total_dat
  ! file name
  character(len=MAX_STRING_LEN),intent(in) :: mesh_file,var_name

  ! local parameters
  integer :: i,j,ie,ier
  real :: val

  ! opens unstructured grid file
  open(IOUT_VTK,file=mesh_file(1:len_trim(mesh_file)),status='unknown',iostat=ier)
  if (ier /= 0 ) then
    print *, 'Error opening VTK output file: ',trim(mesh_file)
    stop 'Error opening VTK output file'
  endif
  write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOUT_VTK,'(a)') 'material model VTK file'
  write(IOUT_VTK,'(a)') 'ASCII'
  write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'

  ! points
  ! writes out all points for each element, not just global ones
  write(IOUT_VTK, '(a,i16,a)') 'POINTS ', ne*8, ' float'
  do ie = 1,ne
    do j = 1,8
      i = total_dat_con(j,ie) + 1  ! needs to add +1 which has been removed before input
      write(IOUT_VTK,'(3e18.6)') real(total_dat_xyz(1,i),kind=4),real(total_dat_xyz(2,i),kind=4),real(total_dat_xyz(3,i),kind=4)
    enddo
  enddo
  write(IOUT_VTK,*) ''

  ! cells
  ! note: indices for VTK start at 0
  write(IOUT_VTK,'(a,i12,i12)') "CELLS ",ne,ne*9
  do ie = 1,ne
    write(IOUT_VTK,'(9i12)') 8,(ie-1)*8,(ie-1)*8+1,(ie-1)*8+2,(ie-1)*8+3, &
                               (ie-1)*8+4,(ie-1)*8+5,(ie-1)*8+6,(ie-1)*8+7
  enddo
  write(IOUT_VTK,*) ''

  ! type: hexahedrons
  write(IOUT_VTK,'(a,i12)') "CELL_TYPES ",ne
  write(IOUT_VTK,'(6i12)') (12,ie = 1,ne)
  write(IOUT_VTK,*) ''

  ! data values
  ! writes out gll-data for each element point
  write(IOUT_VTK,'(a,i12)') "POINT_DATA ",ne*8
  write(IOUT_VTK,'(a)') "SCALARS "//trim(var_name)//" float"
  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
  if (CUSTOM_REAL == SIZE_DOUBLE) then
    ! double precision values
    do ie = 1,ne
      do j = 1,8
        ! converts to float
        i = total_dat_con(j,ie) + 1 ! needs to add +1 which has been removed before input
        val = real(total_dat(i),kind=4)
        ! stay within boundaries of float values, otherwise paraview will complain
        if (abs(val) < tiny(val)) val = sign(1.0,val) * tiny(val)
        if (abs(val) > huge(val)) val = sign(1.0,val) * huge(val)
        write(IOUT_VTK,*) val
      enddo
    enddo
  else
    ! single precision
    do ie = 1,ne
      do j = 1,8
        i = total_dat_con(j,ie) + 1 ! needs to add +1 which has been removed before input
        write(IOUT_VTK,*) total_dat(i)
      enddo
    enddo
  endif
  write(IOUT_VTK,*) ''
  close(IOUT_VTK)

  end subroutine write_VTK_movie_data_elemental


!
!-------------------------------------------------------------------------------------------------
!

! Routine below is commented out since convert='BIG_ENDIAN' is non-standard Fortran (though supported by ifort, gfortran,..).
! It will lead to a compilation warning and error with `make test`


!  subroutine write_VTK_movie_data_binary(ne,np,total_dat_xyz,total_dat_con,total_dat,mesh_file,var_name)
!
!! saves vtk file as binary file for CUSTOM_REAL values
!
!  use constants, only: CUSTOM_REAL,MAX_STRING_LEN,IOUT_VTK
!
!  implicit none
!
!  integer,intent(in) :: ne,np
!
!  ! coordinates
!  real(kind=CUSTOM_REAL), dimension(3,np),intent(in) :: total_dat_xyz
!  ! connections
!  integer, dimension(8,ne),intent(in) :: total_dat_con
!  ! data values array
!  real(kind=CUSTOM_REAL), dimension(np),intent(in) :: total_dat
!  ! file name
!  character(len=MAX_STRING_LEN),intent(in) :: mesh_file,var_name
!
!  ! local parameters
!  integer :: i,it,ier
!  character(len=16) :: str1
!  character(len=12) :: str2,str3
!  character(len=1),parameter :: LF = achar(10) ! line feed character
!
!  ! endianness - determine endianness at run time:
!  ! https://www.ibm.com/developerworks/aix/library/au-endianc/index.html
!  !
!  ! for the Fortran version:
!  ! given the integer value of 1, big endian uses 00 00 00 01 and little endian uses 01 00 00 00 as bit representation.
!  ! using transfer(1,'a') takes the bit-representation of integer value 1 (multi-byte, usually 32-bit)
!  ! and interprets it as a character type (of 1-byte length).
!  ! thus, looking at the first byte, it is either 0 or 1, respectively.
!  ! finally, ichar(c) takes a character and returns its position number.
!  logical, parameter :: is_big_endian = (ichar(transfer(1,'a')) == 0)
!
!  ! note: VTK by default assumes binary data is big endian for "legacy" VTK format (not XML format).
!  !       however, using an intel machine would produce little endian files.
!  !       thus, the convert='BIG_ENDIAN' specifier helps when using an ifort or gfortran compiler.
!  !       unfortunately, it seems not to be supported for cray/ibm compilers...
!
!  !debug
!  !print *,'is big endian: ',is_big_endian
!
!  ! numbers as strings
!  write(str1(1:16),'(i16)') np
!  write(str2(1:12),'(i12)') ne
!  write(str3(1:12),'(i12)') ne*9
!
!  ! opens unstructured grid file as binary file
!  if (is_big_endian) then
!    ! big-endian system
!    open(IOUT_VTK,file=mesh_file(1:len_trim(mesh_file)),access='stream',form='unformatted', &
!         status='unknown',action='write',iostat=ier) ! convert='BIG_ENDIAN' not needed
!  else
!    ! little-endian system
!    open(IOUT_VTK,file=mesh_file(1:len_trim(mesh_file)),access='stream',form='unformatted', &
!         status='unknown',action='write',iostat=ier,convert='BIG_ENDIAN') ! convert='BIG_ENDIAN' specifier is non-standard Fortran
!  endif
!  if (ier /= 0 ) then
!    print *, 'Error opening VTK output file: ',trim(mesh_file)
!    stop 'Error opening VTK output file'
!  endif
!
!  ! file header
!  write(IOUT_VTK) '# vtk DataFile Version 3.1' // LF
!  write(IOUT_VTK) 'material model VTK file' // LF
!  write(IOUT_VTK) 'BINARY' // LF
!  write(IOUT_VTK) 'DATASET UNSTRUCTURED_GRID' // LF
!
!  ! points
!  write(IOUT_VTK) 'POINTS '// str1 // ' float' // LF
!  do i = 1,np
!    write(IOUT_VTK) real(total_dat_xyz(1,i),kind=4),real(total_dat_xyz(2,i),kind=4),real(total_dat_xyz(3,i),kind=4)
!  enddo
!  write(IOUT_VTK) LF
!
!  ! cells
!  ! note: indices for VTK start at 0
!  write(IOUT_VTK) 'CELLS '// str2 // str3 // LF
!  do i = 1,ne
!    write(IOUT_VTK) 8,total_dat_con(1,i),total_dat_con(2,i),total_dat_con(3,i),total_dat_con(4,i), &
!                               total_dat_con(5,i),total_dat_con(6,i),total_dat_con(7,i),total_dat_con(8,i)
!  enddo
!  write(IOUT_VTK) LF
!
!  ! type: hexahedrons
!  write(IOUT_VTK) 'CELL_TYPES ' // str2 // LF
!  write(IOUT_VTK) (12,it = 1,ne)
!  write(IOUT_VTK) LF
!
!  ! data values
!  write(IOUT_VTK) 'POINT_DATA ' // str1 // LF
!  write(IOUT_VTK) 'SCALARS ' // trim(var_name) // ' float' // LF
!  write(IOUT_VTK) 'LOOKUP_TABLE default' // LF
!  do i = 1,np
!    write(IOUT_VTK) real(total_dat(i),kind=4)
!  enddo
!  write(IOUT_VTK) LF
!  close(IOUT_VTK)
!
!  end subroutine write_VTK_movie_data_binary
!
!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_VTU_movie_data(ne,np,total_dat_xyz,total_dat_con,total_dat,mesh_file,var_name)

! saves vtu file in ascii format, for CUSTOM_REAL data values

  use constants, only: CUSTOM_REAL,MAX_STRING_LEN,IOUT_VTK

  implicit none

  integer,intent(in) :: ne,np

  ! coordinates
  real(kind=CUSTOM_REAL), dimension(3,np),intent(in) :: total_dat_xyz
  ! connections
  integer, dimension(8,ne),intent(in) :: total_dat_con
  ! data values array
  real(kind=CUSTOM_REAL), dimension(np),intent(in) :: total_dat
  ! file name
  character(len=MAX_STRING_LEN),intent(in) :: mesh_file,var_name

  ! local parameters
  integer :: i,it,ier
  character(len=16) :: str1
  character(len=12) :: str2

  ! numbers as strings
  write(str1,'(i16)') np
  write(str2,'(i12)') ne

  ! opens unstructured grid file as ascii file
  open(IOUT_VTK,file=mesh_file(1:len_trim(mesh_file)),status='unknown',action='write',iostat=ier)
  if (ier /= 0 ) then
    print *, 'Error opening VTU output file: ',trim(mesh_file)
    stop 'Error opening VTU output file'
  endif

  ! XML file header
  ! see: https://www.vtk.org/Wiki/VTK_XML_Formats
  write(IOUT_VTK,'(a)') '<?xml version="1.0"?>'! avoid white space at beginning of line
  write(IOUT_VTK,'(a)') '<VTKFile type="UnstructuredGrid" version="0.1">'
  write(IOUT_VTK,'(a)') '<UnstructuredGrid>'
  write(IOUT_VTK,'(a)') '<Piece NumberOfPoints="'// trim(adjustl(str1)) // '" NumberOfCells="' // trim(adjustl(str2)) // '">'

  ! points
  write(IOUT_VTK,'(a)') '<Points>'
  write(IOUT_VTK,'(a)') '<DataArray type="Float32" Name="Points" NumberOfComponents="3" format="ascii">'
  do i = 1,np
    write(IOUT_VTK,*) real(total_dat_xyz(1,i),kind=4),real(total_dat_xyz(2,i),kind=4),real(total_dat_xyz(3,i),kind=4)
  enddo
  write(IOUT_VTK,'(a)') '</DataArray>'
  write(IOUT_VTK,'(a)') '</Points>'

  ! cells
  write(IOUT_VTK,'(a)') '<Cells>'
  ! connectivity
  write(IOUT_VTK,'(a)') '<DataArray type="Int64" Name="connectivity" format="ascii">'
  ! note: indices for VTK start at 0
  do i = 1,ne
    write(IOUT_VTK,*) total_dat_con(1,i),total_dat_con(2,i),total_dat_con(3,i),total_dat_con(4,i), &
                      total_dat_con(5,i),total_dat_con(6,i),total_dat_con(7,i),total_dat_con(8,i)
  enddo
  write(IOUT_VTK,'(a)') '</DataArray>'
  ! offsets
  write(IOUT_VTK,'(a)') '<DataArray type="Int64" Name="offsets" format="ascii">'
  write(IOUT_VTK,*) ((it*8),it = 1,ne)
  write(IOUT_VTK,'(a)') '</DataArray>'

  ! type: hexahedrons
  write(IOUT_VTK,'(a)') '<DataArray type="UInt8" Name="types" format="ascii">'
  write(IOUT_VTK,*) (12,it = 1,ne)
  write(IOUT_VTK,'(a)') '</DataArray>'
  write(IOUT_VTK,'(a)') '</Cells>'

  ! empty cell data
  write(IOUT_VTK,'(a)') '<CellData>' // '</CellData>'

  ! data values
  write(IOUT_VTK,'(a)') '<PointData Scalars="Scalars_">'
  write(IOUT_VTK,'(a)') '<DataArray type="Float32" Name="' // trim(var_name) // '" format="ascii">'
  do i = 1,np
    write(IOUT_VTK,*) real(total_dat(i),kind=4)
  enddo
  write(IOUT_VTK,'(a)') '</DataArray>'
  write(IOUT_VTK,'(a)') '</PointData>'

  ! finishes XML file
  write(IOUT_VTK,'(a)') '</Piece>'
  write(IOUT_VTK,'(a)') '</UnstructuredGrid>'

  ! in case of appended data format, with offsets:
  !write(IOUT_VTK,'(a)') '<AppendedData encoding="raw">'
  !write(IOUT_VTK,'(a)') '_'
  !write(IOUT_VTK,*) (real(total_dat_xyz(1,i),kind=4),real(total_dat_xyz(2,i),kind=4),real(total_dat_xyz(3,i),kind=4),i=1,np)
  !write(IOUT_VTK,'(a)') '</AppendedData>'

  write(IOUT_VTK,'(a)') '</VTKFile>'
  close(IOUT_VTK)

  end subroutine write_VTU_movie_data


!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_VTU_movie_data_binary(ne,np,total_dat_xyz,total_dat_con,total_dat,mesh_file,var_name)

! saves vtu file in binary format, for CUSTOM_REAL values

  use constants, only: CUSTOM_REAL,MAX_STRING_LEN,IOUT_VTK,SIZE_REAL,SIZE_INTEGER

  implicit none

  integer,intent(in) :: ne,np

  ! coordinates
  real(kind=CUSTOM_REAL), dimension(3,np),intent(in) :: total_dat_xyz
  ! connections
  integer, dimension(8,ne),intent(in) :: total_dat_con
  ! data values array
  real(kind=CUSTOM_REAL), dimension(np),intent(in) :: total_dat
  ! file name
  character(len=MAX_STRING_LEN),intent(in) :: mesh_file,var_name

  ! local parameters
  integer :: i,it,ier
  integer :: len_bytes,offset
  character(len=16) :: str1,str_endian
  character(len=12) :: str2,str3
  character(len=1),parameter :: LF = achar(10) ! line feed character

  ! endianness - determine endianness at run time:
  ! https://www.ibm.com/developerworks/aix/library/au-endianc/index.html
  !
  ! for the Fortran version:
  ! given the integer value of 1, big endian uses 00 00 00 01 and little endian uses 01 00 00 00 as bit representation.
  ! using transfer(1,'a') takes the bit-representation of integer value 1 (multi-byte, usually 32-bit)
  ! and interprets it as a character type (of 1-byte length).
  ! thus, looking at the first byte, it is either 0 or 1, respectively.
  ! finally, ichar(c) takes a character and returns its position number.
  logical, parameter :: is_big_endian = (ichar(transfer(1,'a')) == 0)

  ! note: VTK by default assumes binary data is big endian for "legacy" VTK format,
  !       we use here the new .vtu file with XML format. in this case, we can specify the endianness by byte_order=".."

  if (is_big_endian) then
    str_endian = 'BigEndian'
  else
    str_endian = 'LittleEndian'
  endif

  ! numbers as strings
  write(str1,'(i16)') np
  write(str2,'(i12)') ne

  ! data offset for appended datasets
  offset = 0
  write(str3,'(i12)') offset

  ! opens unstructured grid file as binary file
  ! convert='BIG_ENDIAN' not needed, will be specified in XML format
  open(IOUT_VTK,file=mesh_file(1:len_trim(mesh_file)),access='stream',form='unformatted', &
         status='unknown',action='write',iostat=ier)
  if (ier /= 0 ) then
    print *, 'Error opening VTU output file: ',trim(mesh_file)
    stop 'Error opening VTU output file'
  endif

  ! XML file header
  ! see: https://www.vtk.org/Wiki/VTK_XML_Formats
  write(IOUT_VTK) '<?xml version="1.0"?>' // LF
  write(IOUT_VTK) '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="'// trim(str_endian) // '">' // LF
  write(IOUT_VTK) '<UnstructuredGrid>' // LF
  write(IOUT_VTK) '<Piece NumberOfPoints="'// str1 // '" NumberOfCells="' // str2 // '">' // LF

  ! points
  write(IOUT_VTK) '<Points>' // LF
  ! binary format: not working properly yet - using appended instead
  !write(IOUT_VTK) '<DataArray type="Float32" Name="Points" NumberOfComponents="3" format="binary" encoding="raw">' // LF
  !! number of bytes to follow
  !! see format description: https://www.vtk.org/Wiki/VTK_XML_Formats#Uncompressed_Data
  !len_bytes = 3 * np * SIZE_REAL
  !write(IOUT_VTK) len_bytes
  !do i = 1,np
  !  write(IOUT_VTK) real(total_dat_xyz(1,i),kind=4),real(total_dat_xyz(2,i),kind=4),real(total_dat_xyz(3,i),kind=4)
  !enddo
  !write(IOUT_VTK) '</DataArray>' // LF
  !
  ! appended format:
  write(IOUT_VTK) '<DataArray type="Float32" Name="Points" NumberOfComponents="3" format="appended" offset="' &
                   // str3 // '"/>' // LF
  ! updates offset
  ! array length in bytes
  len_bytes = 3 * np * SIZE_REAL
  ! new offset: data array size (len_bytes) + 4 bytes for length specifier at the beginning
  offset = offset + len_bytes + 4
  write(str3,'(i12)') offset
  write(IOUT_VTK) '</Points>'//LF

  ! cells
  write(IOUT_VTK) '<Cells>' // LF
  ! connectivity
  write(IOUT_VTK) '<DataArray type="Int32" Name="connectivity" format="appended" offset="' // str3 // '"/>' // LF
  ! updates offset
  len_bytes = 8 * ne * SIZE_INTEGER
  offset = offset + len_bytes + 4
  write(str3,'(i12)') offset

  ! offsets
  write(IOUT_VTK) '<DataArray type="Int32" Name="offsets" format="appended" offset="' // str3 // '"/>' // LF
  ! updates offset
  len_bytes = ne * SIZE_INTEGER
  offset = offset + len_bytes + 4
  write(str3,'(i12)') offset

  ! type: hexahedrons
  write(IOUT_VTK) '<DataArray type="Int32" Name="types" format="appended" offset="' // str3 // '"/>' // LF
  ! updates offset
  len_bytes = ne * SIZE_INTEGER
  offset = offset + len_bytes + 4
  write(str3,'(i12)') offset
  write(IOUT_VTK) '</Cells>' // LF

  ! empty cell data
  write(IOUT_VTK) '<CellData>' // '</CellData>' // LF
  ! data values
  write(IOUT_VTK) '<PointData Scalars="Scalars_">' // LF
  write(IOUT_VTK) '<DataArray type="Float32" Name="' // trim(var_name) // '" format="appended" offset="' &
                  // str3 // '"/>' // LF
  ! updates offset
  !len_bytes = np * SIZE_REAL
  !offset = offset + len_bytes + 4
  !write(str3,'(i12)') offset

  write(IOUT_VTK) '</PointData>' // LF

  ! finishes XML file
  write(IOUT_VTK) '</Piece>' // LF
  write(IOUT_VTK) '</UnstructuredGrid>' // LF

  ! in case of appended data format, with offsets:
  !write(IOUT_VTK) '<AppendedData encoding="base64">' // LF
  write(IOUT_VTK) '<AppendedData encoding="raw">' // LF
  write(IOUT_VTK) '_'
  ! points
  len_bytes = 3 * np * SIZE_REAL
  write(IOUT_VTK) len_bytes
  do i = 1,np
    write(IOUT_VTK) real(total_dat_xyz(1,i),kind=4),real(total_dat_xyz(2,i),kind=4),real(total_dat_xyz(3,i),kind=4)
  enddo
  ! cells
  ! connectivity
  ! number of bytes to follow
  len_bytes = 8 * ne * SIZE_INTEGER
  write(IOUT_VTK) len_bytes
  ! note: indices for VTK start at 0
  do i = 1,ne
    write(IOUT_VTK) total_dat_con(1,i),total_dat_con(2,i),total_dat_con(3,i),total_dat_con(4,i), &
                    total_dat_con(5,i),total_dat_con(6,i),total_dat_con(7,i),total_dat_con(8,i)
  enddo
  ! offsets
  ! number of bytes to follow
  len_bytes = ne * SIZE_INTEGER
  write(IOUT_VTK) len_bytes
  ! offsets (8 points each)
  write(IOUT_VTK) ((it*8),it = 1,ne)
  ! types
  ! number of bytes to follow
  len_bytes = ne * SIZE_INTEGER
  write(IOUT_VTK) len_bytes
  ! type for hex elements is 12
  write(IOUT_VTK) (12,it = 1,ne)

  ! point data values
  ! number of bytes to follow
  len_bytes = np * SIZE_REAL
  write(IOUT_VTK) len_bytes
  ! data values
  do i = 1,np
    write(IOUT_VTK) real(total_dat(i),kind=4)
  enddo

  write(IOUT_VTK) '</AppendedData>' // LF
  write(IOUT_VTK) '</VTKFile>' // LF
  close(IOUT_VTK)

  end subroutine write_VTU_movie_data_binary


!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_VTK_2Ddata_dp(nspec,nglob,xstore_dummy,ystore_dummy,zstore_dummy, &
                                 ibool2D,dat_value,filename)

! stores values on a spherical 2D quad mesh
!
! for details:
!
! VTK type definition numbers:
!   https://vtk.org/doc/nightly/html/vtkCellType_8h_source.html
! ordering images see e.g.:
!   https://lorensen.github.io/VTKExamples/site/VTKBook/05Chapter5/

  use constants, only: IOUT_VTK,MAX_STRING_LEN

  implicit none

  integer :: nspec,nglob

  ! global coordinates
  double precision, dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy
  integer, dimension(4,nspec) :: ibool2D

  ! element flag array
  double precision, dimension(nglob) :: dat_value

  ! file name
  character(len=MAX_STRING_LEN) :: filename

  ! local parameters
  integer :: ispec,i,itype

  ! VTK_QUAD == 8 type, NGNOD2D == 4 corners only
  itype = 8

  ! debug
  !print *,'debug:',itype,nglob,nspec,xstore_dummy(1),ystore_dummy(1),zstore_dummy(1)
  !print *,'debug:',ibool2D(:,1)

  ! write source and receiver VTK files for Paraview
  open(IOUT_VTK,file=trim(filename)//'.vtk',status='unknown')
  write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOUT_VTK,'(a)') 'material model VTK file'
  write(IOUT_VTK,'(a)') 'ASCII'
  write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOUT_VTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do i = 1,nglob
    write(IOUT_VTK,'(3e18.6)') real(xstore_dummy(i),kind=4),real(ystore_dummy(i),kind=4),real(zstore_dummy(i),kind=4)
  enddo
  write(IOUT_VTK,*) ""

  ! note: indices for vtk start at 0
  write(IOUT_VTK,'(a,i12,i12)') "CELLS ",nspec,nspec*(4+1)
  do ispec = 1,nspec
    ! quad4 element using an indexing (left,bottom),(right,bottom),(right,top),(left,top)
    write(IOUT_VTK,'(5i12)') 4,ibool2D(1,ispec)-1,ibool2D(2,ispec)-1,ibool2D(3,ispec)-1,ibool2D(4,ispec)-1
  enddo
  write(IOUT_VTK,*) ""

  ! type: hexahedrons
  write(IOUT_VTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOUT_VTK,'(6i12)') (itype,ispec=1,nspec)
  write(IOUT_VTK,*) ""

  write(IOUT_VTK,'(a,i12)') "POINT_DATA ",nglob
  write(IOUT_VTK,'(a)') "SCALARS dat float"
  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
  do i = 1,nglob
    write(IOUT_VTK,*) real(dat_value(i),kind=4)
  enddo
  write(IOUT_VTK,*) ""
  close(IOUT_VTK)

  end subroutine write_VTK_2Ddata_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_VTU_2Ddata_dp(nspec,nglob,xstore_dummy,ystore_dummy,zstore_dummy, &
                                 ibool2D,dat_value,filename)

! saves vtu file in binary format, for double precision values
!
! stores values on a spherical 2D quad mesh
!
! for details:
!
! VTK type definition numbers:
!   https://vtk.org/doc/nightly/html/vtkCellType_8h_source.html
! ordering images see e.g.:
!   https://lorensen.github.io/VTKExamples/site/VTKBook/05Chapter5/

  use constants, only: IOUT_VTK,MAX_STRING_LEN,SIZE_REAL,SIZE_INTEGER

  implicit none

  integer :: nspec,nglob

  ! global coordinates
  double precision, dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy
  integer, dimension(4,nspec) :: ibool2D

  ! element flag array
  double precision, dimension(nglob) :: dat_value

  ! file name
  character(len=MAX_STRING_LEN) :: filename

  ! local parameters
  integer :: ispec,i,itype,ier
  integer :: len_bytes,offset
  character(len=16) :: str1,str_endian,var_name
  character(len=16) :: str2,str3
  character(len=1),parameter :: LF = achar(10) ! line feed character

  ! endianness - determine endianness at run time:
  logical, parameter :: is_big_endian = (ichar(transfer(1,'a')) == 0)

  ! note: VTK by default assumes binary data is big endian for "legacy" VTK format,
  !       we use here the new .vtu file with XML format. in this case, we can specify the endianness by byte_order=".."
  if (is_big_endian) then
    str_endian = 'BigEndian'
  else
    str_endian = 'LittleEndian'
  endif

  ! VTK_QUAD == 8 type, NGNOD2D == 4 corners only
  itype = 8

  ! array name
  var_name = "val"

  ! integer 32-bit assumed below
  if (SIZE_INTEGER /= 4) stop 'Error vtu output for long int not implemented yet'

  ! integer 32-bit limitation (integer overflow)
  if (nspec > int(2147483646.d0 / dble(4)) .or. nglob > 2147483646.d0 / (3 * dble(SIZE_REAL)) ) then
    print *,'Error: integer 32-bit overflow: nspec = ',nspec,' limit at ',int(2147483646.d0 / dble(4))
    print *,'                                nglob = ',nglob,' limit at ',int(2147483646.d0 / (3*dble(4)))
    stop 'Error vtu data might exceed integer limit'
  endif

  ! numbers as strings
  write(str1,'(i16)') nglob
  write(str2,'(i16)') nspec

  ! data offset for appended datasets
  offset = 0
  write(str3,'(i16)') offset

  ! opens unstructured grid file as binary file
  open(IOUT_VTK,file=trim(filename)//'.vtu',access='stream',form='unformatted', &
       status='unknown',action='write',iostat=ier)
  if (ier /= 0 ) then
    print *, 'Error opening VTU output file: ',trim(filename)
    stop 'Error opening VTU output file'
  endif

  ! XML file header
  ! see: https://www.vtk.org/Wiki/VTK_XML_Formats
  write(IOUT_VTK) '<?xml version="1.0"?>' // LF
  write(IOUT_VTK) '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="'// trim(str_endian) // '">' // LF
  write(IOUT_VTK) '<UnstructuredGrid>' // LF
  write(IOUT_VTK) '<Piece NumberOfPoints="'// trim(adjustl(str1)) // '" NumberOfCells="' // trim(adjustl(str2)) // '">' // LF

  ! points
  write(IOUT_VTK) '<Points>' // LF
  ! binary format: not working properly yet - using appended instead
  ! appended format:
  write(IOUT_VTK) '<DataArray type="Float32" Name="Points" NumberOfComponents="3" format="appended" offset="' &
                   // trim(adjustl(str3)) // '"/>' // LF
  ! updates offset
  ! array length in bytes
  len_bytes = 3 * nglob * SIZE_REAL
  ! new offset: data array size (len_bytes) + 4 bytes for length specifier at the beginning
  offset = offset + len_bytes + 4
  write(str3,'(i16)') offset
  write(IOUT_VTK) '</Points>'//LF

  ! cells
  write(IOUT_VTK) '<Cells>' // LF
  ! connectivity
  write(IOUT_VTK) '<DataArray type="Int32" Name="connectivity" format="appended" offset="' // trim(adjustl(str3)) // '"/>' // LF
  ! updates offset
  len_bytes = 4 * nspec * SIZE_INTEGER
  offset = offset + len_bytes + 4
  write(str3,'(i16)') offset

  ! offsets
  write(IOUT_VTK) '<DataArray type="Int32" Name="offsets" format="appended" offset="' // trim(adjustl(str3)) // '"/>' // LF
  ! updates offset
  len_bytes = nspec * SIZE_INTEGER
  offset = offset + len_bytes + 4
  write(str3,'(i16)') offset

  ! type: hexahedrons
  write(IOUT_VTK) '<DataArray type="Int32" Name="types" format="appended" offset="' // trim(adjustl(str3)) // '"/>' // LF
  ! updates offset
  len_bytes = nspec * SIZE_INTEGER
  offset = offset + len_bytes + 4
  write(str3,'(i16)') offset
  write(IOUT_VTK) '</Cells>' // LF

  ! empty cell data
  write(IOUT_VTK) '<CellData>' // '</CellData>' // LF

  ! point data values
  write(IOUT_VTK) '<PointData Scalars="Scalars_">' // LF
  write(IOUT_VTK) '<DataArray type="Float32" Name="' // trim(var_name) // '" format="appended" offset="' &
                  // trim(adjustl(str3)) // '"/>' // LF
  write(IOUT_VTK) '</PointData>' // LF

  ! finishes XML file
  write(IOUT_VTK) '</Piece>' // LF
  write(IOUT_VTK) '</UnstructuredGrid>' // LF

  ! in case of appended data format, with offsets:
  !write(IOUT_VTK) '<AppendedData encoding="base64">' // LF
  write(IOUT_VTK) '<AppendedData encoding="raw">' // LF
  write(IOUT_VTK) '_'
  ! points
  len_bytes = 3 * nglob * SIZE_REAL
  write(IOUT_VTK) len_bytes
  do i = 1,nglob
    write(IOUT_VTK) real(xstore_dummy(i),kind=4),real(ystore_dummy(i),kind=4),real(zstore_dummy(i),kind=4)
  enddo
  ! cells
  ! number of bytes to follow
  len_bytes = 4 * nspec * SIZE_INTEGER
  write(IOUT_VTK) len_bytes
  ! note: indices for VTK start at 0
  do ispec = 1,nspec
    write(IOUT_VTK) ibool2D(1,ispec)-1,ibool2D(2,ispec)-1,ibool2D(3,ispec)-1,ibool2D(4,ispec)-1
  enddo
  ! offsets
  ! number of bytes to follow
  len_bytes = nspec * SIZE_INTEGER
  write(IOUT_VTK) len_bytes
  ! offsets (4 points each)
  write(IOUT_VTK) (i*4,i = 1,nspec)
  ! types
  ! number of bytes to follow
  len_bytes = nspec * SIZE_INTEGER
  write(IOUT_VTK) len_bytes
  ! type for quads
  write(IOUT_VTK) (itype,i = 1,nspec)
  ! point data values
  ! number of bytes to follow
  len_bytes = nglob * SIZE_REAL
  write(IOUT_VTK) len_bytes
  ! data values
  do i = 1,nglob
    write(IOUT_VTK) real(dat_value(i),kind=4)
  enddo

  write(IOUT_VTK) '</AppendedData>' // LF
  write(IOUT_VTK) '</VTKFile>' // LF
  close(IOUT_VTK)

  end subroutine write_VTU_2Ddata_dp

