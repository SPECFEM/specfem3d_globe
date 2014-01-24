!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            August 2013
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


  program create_slice_VTK

! this programs creates a vtk file that specifies the (velocity) model values on each element point,
! rather than on global points. this will lead to differences especially where a (velocity) discontinuity
! from one element to the other exists. in such cases, this output file should be more accurate
! in how it is visualized.
!
! creates for each slice and each region a new vtk-file, doesn't combine different slices into one file
!
! for compilation, this file 'create_slice_VTK.f90' has to be in the package root directory SPECFEM3D_GLOBE/
!
! cd to your SPECFEM3D_GLOBE root directory:
!   > cd SPECFEM3D_GLOBE/
! create symbolic link:
!   > ln -s UTILS/Visualization/Paraview/create_slice_VTK.f90
! compile with:
!   > f90 -o xcreate_slice_VTK create_slice_VTK.f90
! run :
!   > ./xcreate_slice_VTK my_slices.txt vp ~/DATABASES_MPI/ ~/DATABASES_MPI ~/OUTPUT_FILES/
!
! (or see usage below)
!
  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer,parameter :: MAX_NUM_NODES = 300
  integer  iregion, ir, irs, ire, ires, pfd, efd
  character(len=256) :: sline, arg(7), filename, in_topo_dir, in_file_dir, outdir
  character(len=256) :: prname_topo, prname_file, dimension_file
  character(len=1038) :: command_name
  character(len=256) :: pt_mesh_file1, pt_mesh_file2, mesh_file, em_mesh_file, data_file, topo_file
  integer, dimension(MAX_NUM_NODES) :: node_list, nspec, nglob, npoint, nelement
  integer iproc, num_node, i,j,k,ispec, ios, it, di, dj, dk
  integer np, ne,  njunk
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: data
  real(kind=CUSTOM_REAL),dimension(NGLOB_CRUST_MANTLE) :: xstore, ystore, zstore
  integer ibool(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE)
  integer num_ibool(NGLOB_CRUST_MANTLE)
  logical mask_ibool(NGLOB_CRUST_MANTLE), HIGH_RESOLUTION_MESH
  real x, y, z, dat
  integer numpoin, iglob, n1, n2, n3, n4, n5, n6, n7, n8
  integer iglob1, iglob2, iglob3, iglob4, iglob5, iglob6, iglob7, iglob8


  ! starts here--------------------------------------------------------------------------------------------------
  do i = 1, 6
    call getarg(i,arg(i))
    if (i < 6 .and. trim(arg(i)) == '') then
      print *, ' '
      print *, ' Usage: xcreate_slice_VTK slice_list filename input_topo_dir input_file_dir output_dir [region]'
      print *, ' '
      print *, '   - slice_list:    file containing slice/proc ids '
      print *, '   - filename:    looks for filename.bin must be array of (NGLLX,NGLLY,NGLLZ,nspec) '
      print *, '   - input_topo_dir:    includes "proc***_array_dims.txt '
      print*,  '   - input_file_dir:    includes "proc****filename.bin '
      print *, '   - output_dir:    output mesh files go to here '
      print *, '   if region is not specified, all 3 regions will be collected, otherwise, only collect regions specified'
      print *, ' '
      stop ' Reenter command line options'
    endif
  enddo

  if (NSPEC_CRUST_MANTLE < NSPEC_OUTER_CORE .or. NSPEC_CRUST_MANTLE < NSPEC_INNER_CORE) &
             stop 'This program needs that NSPEC_CRUST_MANTLE > NSPEC_OUTER_CORE and NSPEC_INNER_CORE'

  ! get region id
  if (trim(arg(6)) == '') then
    iregion  = 0
  else
    read(arg(6),*) iregion
  endif
  if (iregion > 3 .or. iregion < 0) stop 'Iregion = 0,1,2,3'
  if (iregion == 0) then
    irs = 1
    ire = 3
  else
    irs = iregion
    ire = irs
  endif

  ! get slices id
  num_node = 0
  open(unit = 20, file = trim(arg(1)), status = 'old',iostat = ios)
  if (ios /= 0) then
    print*,'no file: ',trim(arg(1))
    stop 'Error opening slices file'
  endif

  do while (1 == 1)
    read(20,'(a)',iostat=ios) sline
    if (ios /= 0) exit
    read(sline,*,iostat=ios) njunk
    if (ios /= 0) exit
    num_node = num_node + 1
    node_list(num_node) = njunk
  enddo
  close(20)
  print *, 'slice list: '
  print *, node_list(1:num_node)
  print *, ' '

  ! file to collect
  filename = arg(2)

  ! input and output dir
  in_topo_dir= arg(3)
  in_file_dir= arg(4)
  outdir = arg(5)


  do ir = irs, ire
    print *, '----------- Region ', ir, '----------------'


    ! figure out total number of points and elements for high-res mesh

    do it = 1, num_node

      iproc = node_list(it)

      print *, 'Reading slice ', iproc
      write(prname_topo,'(a,i6.6,a,i1,a)') trim(in_topo_dir)//'/proc',iproc,'_reg',ir,'_'
      dimension_file = trim(prname_topo) //'array_dims.txt'
      open(unit = 27,file = trim(dimension_file),status='old',action='read', iostat = ios)
      if (ios /= 0) then
       print*,'error ',ios
       print*,'file:',trim(dimension_file)
       stop 'Error opening file'
      endif
      read(27,*) nspec(it)
      read(27,*) nglob(it)
      close(27)


    enddo

    print *, 'nspec(it) = ', nspec(1:num_node)
    print *, 'nglob(it) = ', nglob(1:num_node)

    np = 0
    ne = 0

    ! write points information
    do it = 1, num_node

      iproc = node_list(it)


      print *, ' '
      print *, 'Reading slice ', iproc
      write(prname_topo,'(a,i6.6,a,i1,a)') trim(in_topo_dir)//'/proc',iproc,'_reg',ir,'_'
      write(prname_file,'(a,i6.6,a,i1,a)') trim(in_file_dir)//'/proc',iproc,'_reg',ir,'_'

      ! filename.bin
      data_file = trim(prname_file) // trim(filename) // '.bin'
      open(unit = 27,file = trim(data_file),status='old',action='read', iostat = ios,form ='unformatted')
      if (ios /= 0) then
       print*,'error ',ios
       print*,'file:',trim(data_file)
       stop 'Error opening file'
      endif

      data(:,:,:,:) = -1.e10
      read(27) data(:,:,:,1:nspec(it))
      close(27)

      print *,trim(data_file)
      print *,'  min/max value: ',minval(data(:,:,:,1:nspec(it))),maxval(data(:,:,:,1:nspec(it)))
      print *

      ! topology file
      topo_file = trim(prname_topo) // 'solver_data_2' // '.bin'
      open(unit = 28,file = trim(topo_file),status='old',action='read', iostat = ios, form='unformatted')
      if (ios /= 0) then
       print*,'error ',ios
       print*,'file:',trim(topo_file)
       stop 'Error opening file'
      endif
      xstore(:) = 0.0
      ystore(:) = 0.0
      zstore(:) = 0.0
      ibool(:,:,:,:) = -1
      read(28) xstore(1:nglob(it))
      read(28) ystore(1:nglob(it))
      read(28) zstore(1:nglob(it))
      read(28) ibool(:,:,:,1:nspec(it))
      close(28)


      write(mesh_file,'(a,i1,a)') trim(outdir)//'/' // 'reg_',ir,'_'//trim(filename)
      print *, trim(mesh_file)

      ! writes out vtk file
      call write_VTK_data_gll_cr(nspec(it),nglob(it), &
              xstore(1:nglob(it)),ystore(1:nglob(it)),zstore(1:nglob(it)),&
              ibool(:,:,:,1:nspec(it)), &
              data(:,:,:,1:nspec(it)),mesh_file)

    enddo  ! all slices for points

  enddo

  print *, 'Done writing slice files'
  print *, ' '


  end program create_slice_VTK

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_VTK_data_gll_cr(nspec,nglob, &
            xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
            gll_data,prname_file)

! external mesh routine for saving vtk files for custom_real values on all gll points

  implicit none

  include "constants.h"

  integer :: nspec,nglob

! global coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy

! gll data values array
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: gll_data

! file name
  character(len=256) prname_file

  integer :: ispec,i,ier

! write source and receiver VTK files for Paraview
  write(IMAIN,*) '  vtk file: '
  write(IMAIN,*) '    ',prname_file(1:len_trim(prname_file))//'.vtk'

  open(IOUT_VTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown')
  write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOUT_VTK,'(a)') 'material model VTK file'
  write(IOUT_VTK,'(a)') 'ASCII'
  write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'

  ! writes out all points for each element, not just global ones
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
  write(IOUT_VTK,*) ""

  ! note: indices for vtk start at 0
  write(IOUT_VTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
  do ispec=1,nspec
    write(IOUT_VTK,'(9i12)') 8,(ispec-1)*8,(ispec-1)*8+1,(ispec-1)*8+2,(ispec-1)*8+3,&
          (ispec-1)*8+4,(ispec-1)*8+5,(ispec-1)*8+6,(ispec-1)*8+7
  enddo
  write(IOUT_VTK,*) ""

  ! type: hexahedrons
  write(IOUT_VTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOUT_VTK,'(6i12)') (12,ispec=1,nspec)
  write(IOUT_VTK,*) ""

  ! writes out gll-data (velocity) for each element point
  write(IOUT_VTK,'(a,i12)') "POINT_DATA ",nspec*8
  write(IOUT_VTK,'(a)') "SCALARS gll_data float"
  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
  do ispec = 1,nspec
    i = ibool(1,1,1,ispec)
    write(IOUT_VTK,'(3e18.6)') gll_data(1,1,1,ispec)

    i = ibool(NGLLX,1,1,ispec)
    write(IOUT_VTK,'(3e18.6)') gll_data(NGLLX,1,1,ispec)

    i = ibool(NGLLX,NGLLY,1,ispec)
    write(IOUT_VTK,'(3e18.6)') gll_data(NGLLX,NGLLY,1,ispec)

    i = ibool(1,NGLLY,1,ispec)
    write(IOUT_VTK,'(3e18.6)') gll_data(1,NGLLY,1,ispec)

    i = ibool(1,1,NGLLZ,ispec)
    write(IOUT_VTK,'(3e18.6)') gll_data(1,1,NGLLZ,ispec)

    i = ibool(NGLLX,1,NGLLZ,ispec)
    write(IOUT_VTK,'(3e18.6)') gll_data(NGLLX,1,NGLLZ,ispec)

    i = ibool(NGLLX,NGLLY,NGLLZ,ispec)
    write(IOUT_VTK,'(3e18.6)') gll_data(NGLLX,NGLLY,NGLLZ,ispec)

    i = ibool(1,NGLLY,NGLLZ,ispec)-1
    write(IOUT_VTK,'(3e18.6)') gll_data(1,NGLLY,NGLLZ,ispec)
  enddo
  write(IOUT_VTK,*) ""

  close(IOUT_VTK)


  end subroutine write_VTK_data_gll_cr

