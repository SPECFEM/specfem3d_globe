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

program combine_vol_data_vtk

  ! outputs vtk-files (ascii format)

  ! combines the database files on several slices.
  ! the local database file needs to have been collected onto the frontend (copy_local_database.pl)

  implicit none

  include 'constants.h'
  include 'OUTPUT_FILES/values_from_mesher.h'

  integer,parameter :: MAX_NUM_NODES = 2000
  integer  iregion, ir, irs, ire, ires
  character(len=256) :: sline, arg(7), filename, in_topo_dir, in_file_dir, outdir
  character(len=256) :: prname_topo, prname_file, dimension_file
  character(len=256) :: mesh_file
  character(len=256) :: data_file, topo_file
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
  integer ier
  ! instead of taking the first value which appears for a global point, average the values
  ! if there are more than one gll points for a global point (points on element corners, edges, faces)
  logical,parameter:: AVERAGE_GLOBALPOINTS = .false.
  integer:: ibool_count(NGLOB_CRUST_MANTLE)
  real(kind=CUSTOM_REAL):: ibool_dat(NGLOB_CRUST_MANTLE)

  ! note:
  !  if one wants to remove the topography and ellipticity distortion, you would run the mesher again
  !  but turning the flags: TOPOGRAPHY and ELLIPTICITY to .false.
  !  then, use those as topo files: proc***_solver_data.bin
  !  of course, this would also work by just turning ELLIPTICITY to .false. so that the CORRECT_ELLIPTICITY below
  !  becomes unneccessary
  !
  ! puts point locations back into a perfectly spherical shape by removing the ellipticity factor;
  ! useful for plotting spherical cuts at certain depths
  logical,parameter:: CORRECT_ELLIPTICITY = .false.
  integer :: nspl
  double precision :: rspl(NR),espl(NR),espl2(NR)
  logical,parameter :: ONE_CRUST = .false. ! if you want to correct a model with one layer only in PREM crust

  integer, dimension(NSPEC_INNER_CORE) :: idoubling_inner_core ! to get rid of fictitious elements in central cube

  ! global point data
  real,dimension(:),allocatable :: total_dat
  real,dimension(:,:),allocatable :: total_dat_xyz
  integer,dimension(:,:),allocatable :: total_dat_con

  ! starts here--------------------------------------------------------------------------------------------------
  do i = 1, 7
    call getarg(i,arg(i))
    if (i < 7 .and. trim(arg(i)) == '') then
      print *, ' '
      print *, ' Usage: xcombine_vol_data slice_list filename input_topo_dir input_file_dir '
      print *, '        output_dir high/low-resolution [region]'
      print *, ' ***** Notice: now allow different input dir for topo and kernel files ******** '
      print *, '   expect to have the topology and filename.bin(NGLLX,NGLLY,NGLLZ,nspec) '
      print *, '   already collected to input_topo_dir and input_file_dir'
      print *, '   output mesh files (filename_points.mesh, filename_elements.mesh) go to output_dir '
      print *, '   give 0 for low resolution and 1 for high resolution'
      print *, '   if region is not specified, all 3 regions will be collected, otherwise, only collect regions specified'
      stop ' Reenter command line options'
    endif
  enddo

  if (NSPEC_CRUST_MANTLE < NSPEC_OUTER_CORE .or. NSPEC_CRUST_MANTLE < NSPEC_INNER_CORE) &
             stop 'This program needs that NSPEC_CRUST_MANTLE > NSPEC_OUTER_CORE and NSPEC_INNER_CORE'

  ! get region id
  if (trim(arg(7)) == '') then
    iregion  = 0
  else
    read(arg(7),*) iregion
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

  ! resolution
  read(arg(6),*) ires
  di = 0; dj = 0; dk = 0
  if (ires == 0) then
    HIGH_RESOLUTION_MESH = .false.
    di = NGLLX-1; dj = NGLLY-1; dk = NGLLZ-1
  else if( ires == 1 ) then
    HIGH_RESOLUTION_MESH = .true.
    di = 1; dj = 1; dk = 1
  else if( ires == 2 ) then
    HIGH_RESOLUTION_MESH = .false.
    di = (NGLLX-1)/2.0; dj = (NGLLY-1)/2.0; dk = (NGLLZ-1)/2.0
  endif
  if( HIGH_RESOLUTION_MESH ) then
    print *, ' high resolution ', HIGH_RESOLUTION_MESH
  else
    print *, ' low resolution ', HIGH_RESOLUTION_MESH
  endif

  ! sets up ellipticity splines in order to remove ellipticity from point coordinates
  if( CORRECT_ELLIPTICITY ) call make_ellipticity(nspl,rspl,espl,espl2,ONE_CRUST)


  do ir = irs, ire
    print *, '----------- Region ', ir, '----------------'

    ! figure out total number of points and elements for high-res mesh

    do it = 1, num_node

      iproc = node_list(it)

      print *, 'Reading slice ', iproc
      write(prname_topo,'(a,i6.6,a,i1,a)') trim(in_topo_dir)//'/proc',iproc,'_reg',ir,'_'
      write(prname_file,'(a,i6.6,a,i1,a)') trim(in_file_dir)//'/proc',iproc,'_reg',ir,'_'


      dimension_file = trim(prname_topo) //'solver_data.bin'
      open(unit = 27,file = trim(dimension_file),status='old',action='read', iostat = ios, form='unformatted')
      if (ios /= 0) then
       print*,'error ',ios
       print*,'file:',trim(dimension_file)
       stop 'Error opening file'
      endif
      read(27) nspec(it)
      read(27) nglob(it)
      close(27)

      ! check
      if( nspec(it) > NSPEC_CRUST_MANTLE ) stop 'error file nspec too big, please check compilation'
      if( nglob(it) > NGLOB_CRUST_MANTLE ) stop 'error file nglob too big, please check compilation'

      if (HIGH_RESOLUTION_MESH) then
        npoint(it) = nglob(it)
        nelement(it) = nspec(it) * (NGLLX-1) * (NGLLY-1) * (NGLLZ-1)
      else if( ires == 0 ) then
        npoint(it) = nglob(it)
        nelement(it) = nspec(it)
      else if (ires == 2 ) then
        npoint(it) = nglob(it)
        nelement(it) = nspec(it) * (NGLLX-1) * (NGLLY-1) * (NGLLZ-1) / 8
      endif

    enddo

    print *, 'nspec(it) = ', nspec(1:num_node)
    print *, 'nglob(it) = ', nglob(1:num_node)

    !call write_integer_fd(efd,sum(nelement(1:num_node)))

    ! VTK
    print *
    print *,'vtk inital total points: ',sum(npoint(1:num_node))
    print *,'vkt inital total elements: ',sum(nelement(1:num_node))
    print *

    ! creates array to hold point data
    allocate(total_dat(sum(npoint(1:num_node))),stat=ier)
    if( ier /= 0 ) stop 'error allocating total_dat array'
    total_dat(:) = 0.0
    allocate(total_dat_xyz(3,sum(npoint(1:num_node))),stat=ier)
    if( ier /= 0 ) stop 'error allocating total_dat_xyz array'
    total_dat_xyz(:,:) = 0.0
    allocate(total_dat_con(8,sum(nelement(1:num_node))),stat=ier)
    if( ier /= 0 ) stop 'error allocating total_dat_con array'
    total_dat_con(:,:) = 0

    np = 0
    ne = 0

    ! write points information
    do it = 1, num_node

      iproc = node_list(it)

      data(:,:,:,:) = -1.e10

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
      read(27,iostat=ios) data(:,:,:,1:nspec(it))
      if( ios /= 0 ) then
        print*,'read error ',ios
        print*,'file:',trim(data_file)
        stop 'error reading data'
      endif
      close(27)

      print *,trim(data_file)
      print *,'  min/max value: ',minval(data(:,:,:,1:nspec(it))),maxval(data(:,:,:,1:nspec(it)))
      print *

      ! topology file
      topo_file = trim(prname_topo) // 'solver_data.bin'
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
      read(28) nspec(it)
      read(28) nglob(it)
      read(28) xstore(1:nglob(it))
      read(28) ystore(1:nglob(it))
      read(28) zstore(1:nglob(it))
      read(28) ibool(:,:,:,1:nspec(it))
      if (ir==3) read(28) idoubling_inner_core(1:nspec(it)) ! flag that can indicate fictitious elements
      close(28)

      print *, trim(topo_file)


      !average data on global points
      ibool_count(:) = 0
      ibool_dat(:) = 0.0
      if( AVERAGE_GLOBALPOINTS ) then
        do ispec=1,nspec(it)
          ! checks if element counts
          if (ir==3 ) then
            ! inner core
            ! nothing to do for fictitious elements in central cube
            if( idoubling_inner_core(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle
          endif
          ! counts and sums global point data
          do k = 1, NGLLZ, dk
            do j = 1, NGLLY, dj
              do i = 1, NGLLX, di
                iglob = ibool(i,j,k,ispec)

                dat = data(i,j,k,ispec)

                ibool_dat(iglob) = ibool_dat(iglob) + dat
                ibool_count(iglob) = ibool_count(iglob) + 1
              enddo
            enddo
          enddo
        enddo
        do iglob=1,nglob(it)
          if( ibool_count(iglob) > 0 ) then
            ibool_dat(iglob) = ibool_dat(iglob)/ibool_count(iglob)
          endif
        enddo
      endif

      mask_ibool(:) = .false.
      num_ibool(:) = 0
      numpoin = 0


      ! write point file
      do ispec=1,nspec(it)
        ! checks if element counts
        if (ir==3 ) then
          ! inner core
          ! nothing to do for fictitious elements in central cube
          if( idoubling_inner_core(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle
        endif

        ! writes out global point data
        do k = 1, NGLLZ, dk
          do j = 1, NGLLY, dj
            do i = 1, NGLLX, di
              iglob = ibool(i,j,k,ispec)
              if( iglob == -1 ) cycle

              ! takes the averaged data value for mesh
              if( AVERAGE_GLOBALPOINTS ) then
                if(.not. mask_ibool(iglob)) then
                  numpoin = numpoin + 1
                  x = xstore(iglob)
                  y = ystore(iglob)
                  z = zstore(iglob)

                  ! remove ellipticity
                  if( CORRECT_ELLIPTICITY ) call reverse_ellipticity(x,y,z,nspl,rspl,espl,espl2)

                  !dat = data(i,j,k,ispec)
                  dat = ibool_dat(iglob)

                  !call write_real_fd(pfd,x)
                  !call write_real_fd(pfd,y)
                  !call write_real_fd(pfd,z)
                  !call write_real_fd(pfd,dat)

                  ! VTK
                  total_dat(np+numpoin) = dat
                  total_dat_xyz(1,np+numpoin) = x
                  total_dat_xyz(2,np+numpoin) = y
                  total_dat_xyz(3,np+numpoin) = z

                  mask_ibool(iglob) = .true.
                  num_ibool(iglob) = numpoin
                endif
              else
                if(.not. mask_ibool(iglob)) then
                  numpoin = numpoin + 1
                  x = xstore(iglob)
                  y = ystore(iglob)
                  z = zstore(iglob)

                  ! remove ellipticity
                  if( CORRECT_ELLIPTICITY ) call reverse_ellipticity(x,y,z,nspl,rspl,espl,espl2)

                  dat = data(i,j,k,ispec)

                  !call write_real_fd(pfd,x)
                  !call write_real_fd(pfd,y)
                  !call write_real_fd(pfd,z)
                  !call write_real_fd(pfd,dat)

                  ! VTK
                  total_dat(np+numpoin) = dat
                  total_dat_xyz(1,np+numpoin) = x
                  total_dat_xyz(2,np+numpoin) = y
                  total_dat_xyz(3,np+numpoin) = z

                  mask_ibool(iglob) = .true.
                  num_ibool(iglob) = numpoin
                endif
              endif
            enddo ! i
          enddo ! j
        enddo ! k
      enddo !ispec


      ! no way to check the number of points for low-res
      if (HIGH_RESOLUTION_MESH ) then
        if( ir==3 ) then
          npoint(it) = numpoin
        else if( numpoin /= npoint(it)) then
          print*,'region:',ir
          print*,'error number of points:',numpoin,npoint(it)
          stop 'different number of points (high-res)'
        endif
      else if (.not. HIGH_RESOLUTION_MESH) then
        npoint(it) = numpoin
      endif

      ! write elements file
      numpoin = 0
      do ispec = 1, nspec(it)
        ! checks if element counts
        if (ir==3 ) then
          ! inner core
          ! fictitious elements in central cube
          if( idoubling_inner_core(ispec) == IFLAG_IN_FICTITIOUS_CUBE) then
            ! connectivity must be given, otherwise element count would be wrong
            ! maps "fictitious" connectivity, element is all with iglob = 1
            !do k = 1, NGLLZ-1, dk
            !  do j = 1, NGLLY-1, dj
            !    do i = 1, NGLLX-1, di
                  !call write_integer_fd(efd,1)
                  !call write_integer_fd(efd,1)
                  !call write_integer_fd(efd,1)
                  !call write_integer_fd(efd,1)
                  !call write_integer_fd(efd,1)
                  !call write_integer_fd(efd,1)
                  !call write_integer_fd(efd,1)
                  !call write_integer_fd(efd,1)
            !    enddo ! i
            !  enddo ! j
            !enddo ! k
            ! takes next element
            cycle
          endif
        endif

        ! writes out element connectivity
        do k = 1, NGLLZ-1, dk
          do j = 1, NGLLY-1, dj
            do i = 1, NGLLX-1, di

              numpoin = numpoin + 1 ! counts elements

              iglob1 = ibool(i,j,k,ispec)
              iglob2 = ibool(i+di,j,k,ispec)
              iglob3 = ibool(i+di,j+dj,k,ispec)
              iglob4 = ibool(i,j+dj,k,ispec)
              iglob5 = ibool(i,j,k+dk,ispec)
              iglob6 = ibool(i+di,j,k+dk,ispec)
              iglob7 = ibool(i+di,j+dj,k+dk,ispec)
              iglob8 = ibool(i,j+dj,k+dk,ispec)

              n1 = num_ibool(iglob1)+np-1
              n2 = num_ibool(iglob2)+np-1
              n3 = num_ibool(iglob3)+np-1
              n4 = num_ibool(iglob4)+np-1
              n5 = num_ibool(iglob5)+np-1
              n6 = num_ibool(iglob6)+np-1
              n7 = num_ibool(iglob7)+np-1
              n8 = num_ibool(iglob8)+np-1

              !call write_integer_fd(efd,n1)
              !call write_integer_fd(efd,n2)
              !call write_integer_fd(efd,n3)
              !call write_integer_fd(efd,n4)
              !call write_integer_fd(efd,n5)
              !call write_integer_fd(efd,n6)
              !call write_integer_fd(efd,n7)
              !call write_integer_fd(efd,n8)

              ! VTK
              ! note: indices for vtk start at 0
              total_dat_con(1,numpoin + ne) = n1
              total_dat_con(2,numpoin + ne) = n2
              total_dat_con(3,numpoin + ne) = n3
              total_dat_con(4,numpoin + ne) = n4
              total_dat_con(5,numpoin + ne) = n5
              total_dat_con(6,numpoin + ne) = n6
              total_dat_con(7,numpoin + ne) = n7
              total_dat_con(8,numpoin + ne) = n8

            enddo ! i
          enddo ! j
        enddo ! k
      enddo ! ispec

      np = np + npoint(it)
      ne = ne + nelement(it)

    enddo  ! all slices for points

    if (np /= sum(npoint(1:num_node)))  stop 'Error: Number of total points are not consistent'
    if (ne /= sum(nelement(1:num_node))) stop 'Error: Number of total elements are not consistent'

    print *
    print *, 'Total number of points: ', np
    print *, 'Total number of elements: ', ne
    print *

    ! VTK
    ! opens unstructured grid file
    write(mesh_file,'(a,i1,a)') trim(outdir)//'/' // 'reg_',ir,'_'//trim(filename)//'.vtk'
    open(IOVTK,file=mesh_file(1:len_trim(mesh_file)),status='unknown',iostat=ios)
    if( ios /= 0 ) stop 'error opening vtk output file'
    write(IOVTK,'(a)') '# vtk DataFile Version 3.1'
    write(IOVTK,'(a)') 'material model VTK file'
    write(IOVTK,'(a)') 'ASCII'
    write(IOVTK,'(a)') 'DATASET UNSTRUCTURED_GRID'

    ! points
    write(IOVTK, '(a,i16,a)') 'POINTS ', np, ' float'
    do i = 1,np
      write(IOVTK,'(3e18.6)') total_dat_xyz(1,i),total_dat_xyz(2,i),total_dat_xyz(3,i)
    enddo
    write(IOVTK,*) ""

    ! cells
    ! note: indices for vtk start at 0
    write(IOVTK,'(a,i12,i12)') "CELLS ",ne,ne*9
    do i = 1,ne
      write(IOVTK,'(9i12)') 8,total_dat_con(1,i),total_dat_con(2,i),total_dat_con(3,i),total_dat_con(4,i), &
                            total_dat_con(5,i),total_dat_con(6,i),total_dat_con(7,i),total_dat_con(8,i)
    enddo
    write(IOVTK,*) ""

    !call close_file_fd(pfd)
    !call close_file_fd(efd)

    ! add the critical piece: total number of points
    !call open_file_fd(trim(pt_mesh_file2)//char(0),pfd)
    !call write_integer_fd(pfd,np)
    !call close_file_fd(pfd)

    !command_name='cat '//trim(pt_mesh_file2)//' '//trim(pt_mesh_file1)//' '//trim(em_mesh_file)//' > '//trim(mesh_file)
    !print *, ' '
    !print *, 'cat mesh files: '
    !print *, trim(command_name)
    !call system(trim(command_name))

    ! VTK
    ! type: hexahedrons
    write(IOVTK,'(a,i12)') "CELL_TYPES ",ne
    write(IOVTK,*) (12,it=1,ne)
    write(IOVTK,*) ""

    write(IOVTK,'(a,i12)') "POINT_DATA ",np
    write(IOVTK,'(a)') "SCALARS "//trim(filename)//" float"
    write(IOVTK,'(a)') "LOOKUP_TABLE default"
    do i = 1,np
        write(IOVTK,*) total_dat(i)
    enddo
    write(IOVTK,*) ""
    close(IOVTK)

    ! free arrays for this region
    deallocate(total_dat,total_dat_xyz,total_dat_con)


    print *,'written: ',trim(mesh_file)
    print *
  enddo

  print *, 'Done writing mesh files'
  print *, ' '


end program combine_vol_data_vtk

!
! ------------------------------------------------------------------------------------------------
!


  subroutine reverse_ellipticity(x,y,z,nspl,rspl,espl,espl2)

  implicit none

  include "constants.h"

  real(kind=CUSTOM_REAL) :: x,y,z
  integer nspl
  double precision rspl(NR),espl(NR),espl2(NR)
  double precision x1,y1,z1

  double precision ell
  double precision r,theta,phi,factor
  double precision cost,p20

  ! gets spherical coordinates
  x1 = x
  y1 = y
  z1 = z
  call xyz_2_rthetaphi_dble(x1,y1,z1,r,theta,phi)

  cost=dcos(theta)
  p20=0.5d0*(3.0d0*cost*cost-1.0d0)

  ! get ellipticity using spline evaluation
  call spline_evaluation(rspl,espl,espl2,nspl,r,ell)

  factor=ONE-(TWO/3.0d0)*ell*p20

  ! removes ellipticity factor
  x = x / factor
  y = y / factor
  z = z / factor

  end subroutine reverse_ellipticity

!
! ------------------------------------------------------------------------------------------------
!

! copy from make_ellipticity.f90 to avoid compiling issues

  subroutine make_ellipticity(nspl,rspl,espl,espl2,ONE_CRUST)

! creates a spline for the ellipticity profile in PREM
! radius and density are non-dimensional

  implicit none

  include "constants.h"

  integer nspl

  logical ONE_CRUST

! radius of the Earth for gravity calculation
  double precision, parameter :: R_EARTH_ELLIPTICITY = 6371000.d0
! radius of the ocean floor for gravity calculation
  double precision, parameter :: ROCEAN_ELLIPTICITY = 6368000.d0

  double precision rspl(NR),espl(NR),espl2(NR)

  integer i
  double precision ROCEAN,RMIDDLE_CRUST,RMOHO,R80,R220,R400,R600,R670, &
                   R771,RTOPDDOUBLEPRIME,RCMB,RICB
  double precision r_icb,r_cmb,r_topddoubleprime,r_771,r_670,r_600
  double precision r_400,r_220,r_80,r_moho,r_middle_crust,r_ocean,r_0
  double precision r(NR),rho(NR),epsilonval(NR),eta(NR)
  double precision radau(NR),z,k(NR),g_a,bom,exponentval,i_rho,i_radau
  double precision s1(NR),s2(NR),s3(NR)
  double precision yp1,ypn

! PREM
  ROCEAN = 6368000.d0
  RMIDDLE_CRUST = 6356000.d0
  RMOHO = 6346600.d0
  R80  = 6291000.d0
  R220 = 6151000.d0
  R400 = 5971000.d0
  R600 = 5771000.d0
  R670 = 5701000.d0
  R771 = 5600000.d0
  RTOPDDOUBLEPRIME = 3630000.d0
  RCMB = 3480000.d0
  RICB = 1221000.d0

! non-dimensionalize
  r_icb = RICB/R_EARTH_ELLIPTICITY
  r_cmb = RCMB/R_EARTH_ELLIPTICITY
  r_topddoubleprime = RTOPDDOUBLEPRIME/R_EARTH_ELLIPTICITY
  r_771 = R771/R_EARTH_ELLIPTICITY
  r_670 = R670/R_EARTH_ELLIPTICITY
  r_600 = R600/R_EARTH_ELLIPTICITY
  r_400 = R400/R_EARTH_ELLIPTICITY
  r_220 = R220/R_EARTH_ELLIPTICITY
  r_80 = R80/R_EARTH_ELLIPTICITY
  r_moho = RMOHO/R_EARTH_ELLIPTICITY
  r_middle_crust = RMIDDLE_CRUST/R_EARTH_ELLIPTICITY
  r_ocean = ROCEAN_ELLIPTICITY/R_EARTH_ELLIPTICITY
  r_0 = 1.d0

  do i=1,163
    r(i) = r_icb*dble(i-1)/dble(162)
  enddo
  do i=164,323
    r(i) = r_icb+(r_cmb-r_icb)*dble(i-164)/dble(159)
  enddo
  do i=324,336
    r(i) = r_cmb+(r_topddoubleprime-r_cmb)*dble(i-324)/dble(12)
  enddo
  do i=337,517
    r(i) = r_topddoubleprime+(r_771-r_topddoubleprime)*dble(i-337)/dble(180)
  enddo
  do i=518,530
    r(i) = r_771+(r_670-r_771)*dble(i-518)/dble(12)
  enddo
  do i=531,540
    r(i) = r_670+(r_600-r_670)*dble(i-531)/dble(9)
  enddo
  do i=541,565
    r(i) = r_600+(r_400-r_600)*dble(i-541)/dble(24)
  enddo
  do i=566,590
    r(i) = r_400+(r_220-r_400)*dble(i-566)/dble(24)
  enddo
  do i=591,609
    r(i) = r_220+(r_80-r_220)*dble(i-591)/dble(18)
  enddo
  do i=610,619
    r(i) = r_80+(r_moho-r_80)*dble(i-610)/dble(9)
  enddo
  do i=620,626
    r(i) = r_moho+(r_middle_crust-r_moho)*dble(i-620)/dble(6)
  enddo
  do i=627,633
    r(i) = r_middle_crust+(r_ocean-r_middle_crust)*dble(i-627)/dble(6)
  enddo
  do i=634,NR
    r(i) = r_ocean+(r_0-r_ocean)*dble(i-634)/dble(6)
  enddo


! use PREM to get the density profile for ellipticity (fine for other 1D reference models)
  do i=1,NR
    call prem_density(r(i),rho(i),ONE_CRUST,RICB,RCMB,RTOPDDOUBLEPRIME, &
      R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)
    radau(i)=rho(i)*r(i)*r(i)
  enddo

  eta(1)=0.0d0

  k(1)=0.0d0

  do i=2,NR
    call intgrl(i_rho,r,1,i,rho,s1,s2,s3)
    call intgrl(i_radau,r,1,i,radau,s1,s2,s3)
    z=(2.0d0/3.0d0)*i_radau/(i_rho*r(i)*r(i))
    eta(i)=(25.0d0/4.0d0)*((1.0d0-(3.0d0/2.0d0)*z)**2.0d0)-1.0d0
    k(i)=eta(i)/(r(i)**3.0d0)
  enddo

  g_a=4.0D0*i_rho
  bom=TWO_PI/(24.0d0*3600.0d0)
  bom=bom/sqrt(PI*GRAV*RHOAV)
  epsilonval(NR)=15.0d0*(bom**2.0d0)/(24.0d0*i_rho*(eta(NR)+2.0d0))

  do i=1,NR-1
    call intgrl(exponentval,r,i,NR,k,s1,s2,s3)
    epsilonval(i)=epsilonval(NR)*exp(-exponentval)
  enddo

! get ready to spline epsilonval
  nspl=1
  rspl(1)=r(1)
  espl(1)=epsilonval(1)
  do i=2,NR
    if(r(i) /= r(i-1)) then
      nspl=nspl+1
      rspl(nspl)=r(i)
      espl(nspl)=epsilonval(i)
    endif
  enddo

! spline epsilonval
  yp1=0.0d0
  ypn=(5.0d0/2.0d0)*(bom**2)/g_a-2.0d0*epsilonval(NR)
  call spline_construction(rspl,espl,nspl,yp1,ypn,espl2)

  end subroutine make_ellipticity

!
! ------------------------------------------------------------------------------------------------
!

! copy from model_prem.f90 to avoid compiling issues

  subroutine prem_density(x,rho,ONE_CRUST,RICB,RCMB,RTOPDDOUBLEPRIME, &
      R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)

  implicit none

  include "constants.h"

  double precision x,rho,RICB,RCMB,RTOPDDOUBLEPRIME, &
      R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN

  logical ONE_CRUST

  double precision r

  ! compute real physical radius in meters
  r = x * R_EARTH

  ! calculates density according to radius
  if(r <= RICB) then
    rho=13.0885d0-8.8381d0*x*x
  else if(r > RICB .and. r <= RCMB) then
    rho=12.5815d0-1.2638d0*x-3.6426d0*x*x-5.5281d0*x*x*x
  else if(r > RCMB .and. r <= RTOPDDOUBLEPRIME) then
    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
  else if(r > RTOPDDOUBLEPRIME .and. r <= R771) then
    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
  else if(r > R771 .and. r <= R670) then
    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
  else if(r > R670 .and. r <= R600) then
    rho=5.3197d0-1.4836d0*x
  else if(r > R600 .and. r <= R400) then
    rho=11.2494d0-8.0298d0*x
  else if(r > R400 .and. r <= R220) then
    rho=7.1089d0-3.8045d0*x
  else if(r > R220 .and. r <= R80) then
    rho=2.6910d0+0.6924d0*x
  else
    if(r > R80 .and. r <= RMOHO) then
      rho=2.6910d0+0.6924d0*x
    else if(r > RMOHO .and. r <= RMIDDLE_CRUST) then
      if(ONE_CRUST) then
        rho=2.6d0
      else
        rho=2.9d0
      endif
    else if(r > RMIDDLE_CRUST .and. r <= ROCEAN) then
      rho=2.6d0
    else if(r > ROCEAN) then
      rho=2.6d0
    endif
  endif

  rho=rho*1000.0d0/RHOAV

  end subroutine prem_density

!
! ------------------------------------------------------------------------------------------------
!

! copy from intgrl.f90 to avoid compiling issues


 subroutine intgrl(sum,r,nir,ner,f,s1,s2,s3)

! Computes the integral of f[i]*r[i]*r[i] from i=nir to i=ner for
! radii values as in model PREM_an640

  implicit none

! Argument variables
  integer ner,nir
  double precision f(640),r(640),s1(640),s2(640)
  double precision s3(640),sum

! Local variables
  double precision, parameter :: third = 1.0d0/3.0d0
  double precision, parameter :: fifth = 1.0d0/5.0d0
  double precision, parameter :: sixth = 1.0d0/6.0d0

  double precision rji,yprime(640)
  double precision s1l,s2l,s3l

  integer i,j,n,kdis(28)
  integer ndis,nir1



  data kdis/163,323,336,517,530,540,565,590,609,619,626,633,16*0/

  ndis = 12
  n = 640

  call deriv(f,yprime,n,r,ndis,kdis,s1,s2,s3)
  nir1 = nir + 1
  sum = 0.0d0
  do i=nir1,ner
    j = i-1
    rji = r(i) - r(j)
    s1l = s1(j)
    s2l = s2(j)
    s3l = s3(j)
    sum = sum + r(j)*r(j)*rji*(f(j) &
              + rji*(0.5d0*s1l + rji*(third*s2l + rji*0.25d0*s3l))) &
              + 2.0d0*r(j)*rji*rji*(0.5d0*f(j) + rji*(third*s1l + rji*(0.25d0*s2l + rji*fifth*s3l))) &
              + rji*rji*rji*(third*f(j) + rji*(0.25d0*s1l + rji*(fifth*s2l + rji*sixth*s3l)))
  enddo

  end subroutine intgrl

! -------------------------------

  subroutine deriv(y,yprime,n,r,ndis,kdis,s1,s2,s3)

  implicit none

! Argument variables
  integer kdis(28),n,ndis
  double precision r(n),s1(n),s2(n),s3(n)
  double precision y(n),yprime(n)

! Local variables
  integer i,j,j1,j2
  integer k,nd,ndp
  double precision a0,b0,b1
  double precision f(3,1000),h,h2,h2a
  double precision h2b,h3a,ha,s13
  double precision s21,s32,yy(3)

  yy(1) = 0.d0
  yy(2) = 0.d0
  yy(3) = 0.d0

  ndp=ndis+1
  do 3 nd=1,ndp
  if(nd == 1) goto 4
  if(nd == ndp) goto 5
  j1=kdis(nd-1)+1
  j2=kdis(nd)-2
  goto 6
    4 j1=1
  j2=kdis(1)-2
  goto 6
    5 j1=kdis(ndis)+1
  j2=n-2
    6 if((j2+1-j1)>0) goto 11
  j2=j2+2
  yy(1)=(y(j2)-y(j1))/(r(j2)-r(j1))
  s1(j1)=yy(1)
  s1(j2)=yy(1)
  s2(j1)=yy(2)
  s2(j2)=yy(2)
  s3(j1)=yy(3)
  s3(j2)=yy(3)
  goto 3
   11 a0=0.0d0
  if(j1 == 1) goto 7
  h=r(j1+1)-r(j1)
  h2=r(j1+2)-r(j1)
  yy(1)=h*h2*(h2-h)
  h=h*h
  h2=h2*h2
  b0=(y(j1)*(h-h2)+y(j1+1)*h2-y(j1+2)*h)/yy(1)
  goto 8
 7 b0=0.0d0
 8 b1=b0

  if(j2 > 1000) stop 'error in subroutine deriv for j2'

  do i=j1,j2
    h=r(i+1)-r(i)
    yy(1)=y(i+1)-y(i)
    h2=h*h
    ha=h-a0
    h2a=h-2.0d0*a0
    h3a=2.0d0*h-3.0d0*a0
    h2b=h2*b0
    s1(i)=h2/ha
    s2(i)=-ha/(h2a*h2)
    s3(i)=-h*h2a/h3a
    f(1,i)=(yy(1)-h*b0)/(h*ha)
    f(2,i)=(h2b-yy(1)*(2.0d0*h-a0))/(h*h2*h2a)
    f(3,i)=-(h2b-3.0d0*yy(1)*ha)/(h*h3a)
    a0=s3(i)
    b0=f(3,i)
  enddo

  i=j2+1
  h=r(i+1)-r(i)
  yy(1)=y(i+1)-y(i)
  h2=h*h
  ha=h-a0
  h2a=h*ha
  h2b=h2*b0-yy(1)*(2.d0*h-a0)
  s1(i)=h2/ha
  f(1,i)=(yy(1)-h*b0)/h2a
  ha=r(j2)-r(i+1)
  yy(1)=-h*ha*(ha+h)
  ha=ha*ha
  yy(1)=(y(i+1)*(h2-ha)+y(i)*ha-y(j2)*h2)/yy(1)
  s3(i)=(yy(1)*h2a+h2b)/(h*h2*(h-2.0d0*a0))
  s13=s1(i)*s3(i)
  s2(i)=f(1,i)-s13

  do j=j1,j2
    k=i-1
    s32=s3(k)*s2(i)
    s1(i)=f(3,k)-s32
    s21=s2(k)*s1(i)
    s3(k)=f(2,k)-s21
    s13=s1(k)*s3(k)
    s2(k)=f(1,k)-s13
    i=k
  enddo

  s1(i)=b1
  j2=j2+2
  s1(j2)=yy(1)
  s2(j2)=yy(2)
  s3(j2)=yy(3)
 3 continue

  do i=1,n
    yprime(i)=s1(i)
  enddo

  end subroutine deriv

!
! ------------------------------------------------------------------------------------------------
!

! copy from spline_routines.f90 to avoid compiling issues

! compute spline coefficients

  subroutine spline_construction(xpoint,ypoint,npoint,tangent_first_point,tangent_last_point,spline_coefficients)

  implicit none

! tangent to the spline imposed at the first and last points
  double precision, intent(in) :: tangent_first_point,tangent_last_point

! number of input points and coordinates of the input points
  integer, intent(in) :: npoint
  double precision, dimension(npoint), intent(in) :: xpoint,ypoint

! spline coefficients output by the routine
  double precision, dimension(npoint), intent(out) :: spline_coefficients

  integer :: i

  double precision, dimension(:), allocatable :: temporary_array

  allocate(temporary_array(npoint))

  spline_coefficients(1) = - 1.d0 / 2.d0

  temporary_array(1) = (3.d0/(xpoint(2)-xpoint(1)))*((ypoint(2)-ypoint(1))/(xpoint(2)-xpoint(1))-tangent_first_point)

  do i = 2,npoint-1

    spline_coefficients(i) = ((xpoint(i)-xpoint(i-1))/(xpoint(i+1)-xpoint(i-1))-1.d0) &
       / ((xpoint(i)-xpoint(i-1))/(xpoint(i+1)-xpoint(i-1))*spline_coefficients(i-1)+2.d0)

    temporary_array(i) = (6.d0*((ypoint(i+1)-ypoint(i))/(xpoint(i+1)-xpoint(i)) &
       - (ypoint(i)-ypoint(i-1))/(xpoint(i)-xpoint(i-1)))/(xpoint(i+1)-xpoint(i-1)) &
       - (xpoint(i)-xpoint(i-1))/(xpoint(i+1)-xpoint(i-1))*temporary_array(i-1)) &
       / ((xpoint(i)-xpoint(i-1))/(xpoint(i+1)-xpoint(i-1))*spline_coefficients(i-1)+2.d0)

  enddo

  spline_coefficients(npoint) = ((3.d0/(xpoint(npoint)-xpoint(npoint-1))) &
      * (tangent_last_point-(ypoint(npoint)-ypoint(npoint-1))/(xpoint(npoint)-xpoint(npoint-1))) &
      - 1.d0/2.d0*temporary_array(npoint-1))/(1.d0/2.d0*spline_coefficients(npoint-1)+1.d0)

  do i = npoint-1,1,-1
    spline_coefficients(i) = spline_coefficients(i)*spline_coefficients(i+1) + temporary_array(i)
  enddo

  deallocate(temporary_array)

  end subroutine spline_construction

! --------------

! evaluate a spline

  subroutine spline_evaluation(xpoint,ypoint,spline_coefficients,npoint,x_evaluate_spline,y_spline_obtained)

  implicit none

! number of input points and coordinates of the input points
  integer, intent(in) :: npoint
  double precision, dimension(npoint), intent(in) :: xpoint,ypoint

! spline coefficients to use
  double precision, dimension(npoint), intent(in) :: spline_coefficients

! abscissa at which we need to evaluate the value of the spline
  double precision, intent(in):: x_evaluate_spline

! ordinate evaluated by the routine for the spline at this abscissa
  double precision, intent(out):: y_spline_obtained

  integer :: index_loop,index_lower,index_higher

  double precision :: coef1,coef2

! initialize to the whole interval
  index_lower = 1
  index_higher = npoint

! determine the right interval to use, by dichotomy
  do while (index_higher - index_lower > 1)
! compute the middle of the interval
    index_loop = (index_higher + index_lower) / 2
    if(xpoint(index_loop) > x_evaluate_spline) then
      index_higher = index_loop
    else
      index_lower = index_loop
    endif
  enddo

! test that the interval obtained does not have a size of zero
! (this could happen for instance in the case of duplicates in the input list of points)
  if(xpoint(index_higher) == xpoint(index_lower)) stop 'incorrect interval found in spline evaluation'

  coef1 = (xpoint(index_higher) - x_evaluate_spline) / (xpoint(index_higher) - xpoint(index_lower))
  coef2 = (x_evaluate_spline - xpoint(index_lower)) / (xpoint(index_higher) - xpoint(index_lower))

  y_spline_obtained = coef1*ypoint(index_lower) + coef2*ypoint(index_higher) + &
        ((coef1**3 - coef1)*spline_coefficients(index_lower) + &
         (coef2**3 - coef2)*spline_coefficients(index_higher))*((xpoint(index_higher) - xpoint(index_lower))**2)/6.d0

  end subroutine spline_evaluation


!
! ------------------------------------------------------------------------------------------------
!

! copy from rthetaphi_xyz.f90 to avoid compiling issues


  subroutine xyz_2_rthetaphi_dble(x,y,z,r,theta,phi)

! convert x y z to r theta phi, double precision call

  implicit none

  include "constants.h"

  double precision x,y,z,r,theta,phi
  double precision xmesh,ymesh,zmesh

  xmesh = x
  ymesh = y
  zmesh = z

  if(zmesh > -SMALL_VAL_ANGLE .and. zmesh <= ZERO) zmesh = -SMALL_VAL_ANGLE
  if(zmesh < SMALL_VAL_ANGLE .and. zmesh >= ZERO) zmesh = SMALL_VAL_ANGLE

  theta = datan2(dsqrt(xmesh*xmesh+ymesh*ymesh),zmesh)

  if(xmesh > -SMALL_VAL_ANGLE .and. xmesh <= ZERO) xmesh = -SMALL_VAL_ANGLE
  if(xmesh < SMALL_VAL_ANGLE .and. xmesh >= ZERO) xmesh = SMALL_VAL_ANGLE

  phi = datan2(ymesh,xmesh)

  r = dsqrt(xmesh*xmesh + ymesh*ymesh + zmesh*zmesh)

  end subroutine xyz_2_rthetaphi_dble

