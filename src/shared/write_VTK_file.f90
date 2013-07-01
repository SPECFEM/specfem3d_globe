!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
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


  subroutine write_VTK_data_points(nglob, &
                                  xstore_dummy,ystore_dummy,zstore_dummy, &
                                  points_globalindices,num_points_globalindices, &
                                  prname_file)

! external mesh routine for saving vtk files for points locations

  implicit none

  include "constants.h"

  integer :: nglob

  ! global coordinates
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy

  ! gll data values array
  integer :: num_points_globalindices
  integer, dimension(num_points_globalindices) :: points_globalindices

  ! file name
  character(len=150) prname_file

  integer :: i,iglob

  ! write source and receiver VTK files for Paraview
  !debug
  !write(IMAIN,*) '  vtk file: '
  !write(IMAIN,*) '    ',prname_file(1:len_trim(prname_file))//'.vtk'

  open(IOVTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown')
  write(IOVTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOVTK,'(a)') 'material model VTK file'
  write(IOVTK,'(a)') 'ASCII'
  write(IOVTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOVTK, '(a,i12,a)') 'POINTS ', num_points_globalindices, ' float'
  do i=1,num_points_globalindices
    iglob = points_globalindices(i)
    if( iglob <= 0 .or. iglob > nglob ) then
      print*,'error: '//prname_file(1:len_trim(prname_file))//'.vtk'
      print*,'error global index: ',iglob,i
      close(IOVTK)
      stop 'error vtk points file'
    endif

    write(IOVTK,'(3e18.6)') xstore_dummy(iglob),ystore_dummy(iglob),zstore_dummy(iglob)
  enddo
  write(IOVTK,*) ""

  close(IOVTK)

  end subroutine write_VTK_data_points

!
!-------------------------------------------------------------------------------------------------
!


  subroutine write_VTK_glob_points(nglob, &
                                  xstore_dummy,ystore_dummy,zstore_dummy, &
                                  glob_values, &
                                  prname_file)

! external mesh routine for saving vtk files for points locations

  implicit none

  include "constants.h"

  integer :: nglob

  ! global coordinates
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy

  ! gll data values array
  real(kind=CUSTOM_REAL), dimension(nglob) :: glob_values

  ! file name
  character(len=150) prname_file

  integer :: iglob

  ! write source and receiver VTK files for Paraview
  !debug
  !write(IMAIN,*) '  vtk file: '
  !write(IMAIN,*) '    ',prname_file(1:len_trim(prname_file))//'.vtk'

  open(IOVTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown')
  write(IOVTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOVTK,'(a)') 'material model VTK file'
  write(IOVTK,'(a)') 'ASCII'
  write(IOVTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOVTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do iglob=1,nglob
    write(IOVTK,*) xstore_dummy(iglob),ystore_dummy(iglob),zstore_dummy(iglob)
  enddo
  write(IOVTK,*) ""

  ! writes out gll-data (velocity) for each element point
  write(IOVTK,'(a,i12)') "POINT_DATA ",nglob
  write(IOVTK,'(a)') "SCALARS glob_data float"
  write(IOVTK,'(a)') "LOOKUP_TABLE default"
  do iglob=1,nglob
    write(IOVTK,*) glob_values(iglob)
  enddo
  write(IOVTK,*) ""

  close(IOVTK)

  end subroutine write_VTK_glob_points



!
!-------------------------------------------------------------------------------------------------
!


  subroutine write_VTK_data_elem_l(nspec,nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                        elem_flag,prname_file)

! routine for saving vtk file holding logical flag on each spectral element

  implicit none

  include "constants.h"

  integer :: nspec,nglob

  ! global coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy

  ! element flag array
  logical, dimension(nspec) :: elem_flag
  integer :: ispec,i

  ! file name
  character(len=150) prname_file

  ! write source and receiver VTK files for Paraview
  !debug
  !write(IMAIN,*) '  vtk file: '
  !write(IMAIN,*) '    ',prname_file(1:len_trim(prname_file))//'.vtk'

  open(IOVTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown')
  write(IOVTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOVTK,'(a)') 'material model VTK file'
  write(IOVTK,'(a)') 'ASCII'
  write(IOVTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOVTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do i=1,nglob
    write(IOVTK,'(3e18.6)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)
  enddo
  write(IOVTK,*) ""

  ! note: indices for vtk start at 0
  write(IOVTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
  do ispec=1,nspec
    write(IOVTK,'(9i12)') 8,ibool(1,1,1,ispec)-1,ibool(NGLLX,1,1,ispec)-1,ibool(NGLLX,NGLLY,1,ispec)-1,ibool(1,NGLLY,1,ispec)-1,&
          ibool(1,1,NGLLZ,ispec)-1,ibool(NGLLX,1,NGLLZ,ispec)-1,ibool(NGLLX,NGLLY,NGLLZ,ispec)-1,ibool(1,NGLLY,NGLLZ,ispec)-1
  enddo
  write(IOVTK,*) ""

  ! type: hexahedrons
  write(IOVTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOVTK,*) (12,ispec=1,nspec)
  write(IOVTK,*) ""

  write(IOVTK,'(a,i12)') "CELL_DATA ",nspec
  write(IOVTK,'(a)') "SCALARS elem_flag integer"
  write(IOVTK,'(a)') "LOOKUP_TABLE default"
  do ispec = 1,nspec
    if( elem_flag(ispec) .eqv. .true. ) then
      write(IOVTK,*) 1
    else
      write(IOVTK,*) 0
    endif
  enddo
  write(IOVTK,*) ""
  close(IOVTK)


  end subroutine write_VTK_data_elem_l


!
!-------------------------------------------------------------------------------------------------
!


  subroutine write_VTK_data_elem_i(nspec,nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                        elem_flag,prname_file)


! routine for saving vtk file holding integer value on each spectral element

  implicit none

  include "constants.h"

  integer :: nspec,nglob

  ! global coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy

  ! element flag array
  integer, dimension(nspec) :: elem_flag
  integer :: ispec,i

  ! file name
  character(len=150) prname_file

  ! write source and receiver VTK files for Paraview
  !debug
  !write(IMAIN,*) '  vtk file: '
  !write(IMAIN,*) '    ',prname_file(1:len_trim(prname_file))//'.vtk'

  open(IOVTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown')
  write(IOVTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOVTK,'(a)') 'material model VTK file'
  write(IOVTK,'(a)') 'ASCII'
  write(IOVTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOVTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do i=1,nglob
    write(IOVTK,'(3e18.6)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)
  enddo
  write(IOVTK,*) ""

  ! note: indices for vtk start at 0
  write(IOVTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
  do ispec=1,nspec
    write(IOVTK,'(9i12)') 8,ibool(1,1,1,ispec)-1,ibool(NGLLX,1,1,ispec)-1,ibool(NGLLX,NGLLY,1,ispec)-1,ibool(1,NGLLY,1,ispec)-1,&
          ibool(1,1,NGLLZ,ispec)-1,ibool(NGLLX,1,NGLLZ,ispec)-1,ibool(NGLLX,NGLLY,NGLLZ,ispec)-1,ibool(1,NGLLY,NGLLZ,ispec)-1
  enddo
  write(IOVTK,*) ""

  ! type: hexahedrons
  write(IOVTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOVTK,*) (12,ispec=1,nspec)
  write(IOVTK,*) ""

  write(IOVTK,'(a,i12)') "CELL_DATA ",nspec
  write(IOVTK,'(a)') "SCALARS elem_val integer"
  write(IOVTK,'(a)') "LOOKUP_TABLE default"
  do ispec = 1,nspec
    write(IOVTK,*) elem_flag(ispec)
  enddo
  write(IOVTK,*) ""
  close(IOVTK)

  end subroutine write_VTK_data_elem_i

!
!-------------------------------------------------------------------------------------------------
!

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
  write(IOVTK,*) ""

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
  write(IOVTK,*) ""

  ! type: hexahedrons
  write(IOVTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOVTK,*) (12,ispec=1,nspec)
  write(IOVTK,*) ""

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
  write(IOVTK,*) ""

  close(IOVTK)


  end subroutine write_VTK_data_cr

!
!-------------------------------------------------------------------------------------------------
!

! external mesh routine for saving vtk files for custom_real values on global points

  subroutine write_VTK_data_cr_all(myrank,NPROCTOT,idoubling, &
                              nspec,nglob, &
                              xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                              glob_data,prname_file)

! outputs single file for all processes

  implicit none

  include "constants.h"

  include 'mpif.h'
  include "precision.h"

  integer :: myrank,NPROCTOT

  integer ::nspec,nglob

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
    allocate(store_val_x_all(nglob,0:NPROCTOT-1), &
            store_val_y_all(nglob,0:NPROCTOT-1), &
            store_val_z_all(nglob,0:NPROCTOT-1), &
            store_val_ux_all(nglob,0:NPROCTOT-1), &
            store_val_uy_all(nglob,0:NPROCTOT-1), &
            store_val_uz_all(nglob,0:NPROCTOT-1), &
            idoubling_all(nspec,0:NPROCTOT-1), &
            ibool_all(NGLLX,NGLLY,NGLLZ,nspec,0:NPROCTOT-1),stat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error allocating stores')
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
    if( ier /= 0 ) call exit_mpi(myrank,'error allocating dummy stores')
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
    write(IOVTK, '(a,i12,a)') 'POINTS ', nglob*NPROCTOT, ' float'
    do iproc=0, NPROCTOT-1
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
    write(IOVTK,*) ""

    ! defines cell on coarse corner points
    ! note: indices for vtk start at 0
    write(IOVTK,'(a,i12,i12)') "CELLS ",nspec*NPROCTOT,nspec*NPROCTOT*9
    do iproc=0, NPROCTOT-1
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
    write(IOVTK,*) ""

    ! type: hexahedrons
    write(IOVTK,'(a,i12)') "CELL_TYPES ",nspec*NPROCTOT
    write(IOVTK,*) (12,ispec=1,nspec*NPROCTOT)
    write(IOVTK,*) ""

    ! x components
    write(IOVTK,'(a,i12)') "POINT_DATA ",nglob*NPROCTOT
    write(IOVTK,'(a)') "SCALARS x_comp float"
    write(IOVTK,'(a)') "LOOKUP_TABLE default"
    do iproc=0, NPROCTOT-1
      do i = 1,nglob
        write(IOVTK,*) store_val_ux_all(i,iproc)
      enddo
    enddo
    ! y components
    write(IOVTK,'(a)') "SCALARS y_comp float"
    write(IOVTK,'(a)') "LOOKUP_TABLE default"
    do iproc=0, NPROCTOT-1
      do i = 1,nglob
        write(IOVTK,*) store_val_uy_all(i,iproc)
      enddo
    enddo
    ! z components
    write(IOVTK,'(a)') "SCALARS z_comp float"
    write(IOVTK,'(a)') "LOOKUP_TABLE default"
    do iproc=0, NPROCTOT-1
      do i = 1,nglob
        write(IOVTK,*) store_val_uz_all(i,iproc)
      enddo
    enddo
    ! norm
    write(IOVTK,'(a)') "SCALARS norm float"
    write(IOVTK,'(a)') "LOOKUP_TABLE default"
    do iproc=0, NPROCTOT-1
      do i = 1,nglob
        write(IOVTK,*) sqrt( store_val_ux_all(i,iproc)**2 &
                          + store_val_uy_all(i,iproc)**2 &
                          + store_val_uz_all(i,iproc)**2 )
      enddo
    enddo
    write(IOVTK,*) ""

    close(IOVTK)

  endif

  deallocate(store_val_x_all,store_val_y_all,store_val_z_all, &
            store_val_ux_all,store_val_uy_all,store_val_uz_all, &
            ibool_all)

  end subroutine write_VTK_data_cr_all
