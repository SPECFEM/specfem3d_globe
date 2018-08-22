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

program combine_vol_data

  ! outputs VTK files (ASCII format)

  use constants

  ! combines the database files on several slices.
  ! the local database file needs to have been collected onto the frontend (copy_local_database.pl)
#ifdef ADIOS_INPUT
  use adios_read_mod
  use combine_vol_data_adios_mod
#endif

  implicit none

  include "OUTPUT_FILES/values_from_mesher.h"

  integer,parameter :: MAX_NUM_NODES = 2000
  integer, dimension(MAX_NUM_NODES) :: node_list, nspec_list, nglob_list
  integer, dimension(MAX_NUM_NODES) :: npoint, nelement

  ! dp: or use NPROCTOT_VAL to avoid setting max_num_nodes?
  !integer, dimension(NPROCTOT_VAL) :: node_list, nspec_list, nglob_list

  integer :: ir, irs, ire, ires
  integer :: iproc, num_node, i,j,k,ispec, it, di, dj, dk
  integer :: np, ne
  integer :: ier

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: data
  real(kind=CUSTOM_REAL),dimension(NGLOB_CRUST_MANTLE) :: xstore, ystore, zstore
  integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool

  integer :: num_ibool(NGLOB_CRUST_MANTLE)
  logical :: mask_ibool(NGLOB_CRUST_MANTLE)
  logical :: HIGH_RESOLUTION_MESH

  real(kind=CUSTOM_REAL) :: x, y, z
  real :: dat
  integer :: numpoin, iglob, n1, n2, n3, n4, n5, n6, n7, n8
  integer :: iglob1, iglob2, iglob3, iglob4, iglob5, iglob6, iglob7, iglob8

  ! instead of taking the first value which appears for a global point, average the values
  ! if there are more than one GLL points for a global point (points on element corners, edges, faces)
  logical,parameter:: AVERAGE_GLOBALPOINTS = .false.

  integer,dimension(:),allocatable :: ibool_count
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: data_avg

  ! note:
  !  if one wants to remove the topography and ellipticity distortion, you would run the mesher again
  !  but turning the flags: TOPOGRAPHY and ELLIPTICITY to .false.
  !  then, use those as topo files: proc***_solver_data.bin
  !  of course, this would also work by just turning ELLIPTICITY to .false. so that the CORRECT_ELLIPTICITY below
  !  becomes unnecessary
  !
  ! puts point locations back into a perfectly spherical shape by removing the ellipticity factor;
  ! useful for plotting spherical cuts at certain depths
  logical, parameter :: CORRECT_ELLIPTICITY = .false.

  integer :: nspl
  double precision :: rspl(NR),espl(NR),espl2(NR)
  logical,parameter :: ONE_CRUST = .false. ! if you want to correct a model with one layer only in PREM crust

  integer, dimension(NSPEC_INNER_CORE) :: idoubling_inner_core ! to get rid of fictitious elements in central cube

  character(len=MAX_STRING_LEN) :: arg(7)
  character(len=MAX_STRING_LEN) :: filename, outdir
  character(len=MAX_STRING_LEN) :: data_file, var_name

#ifdef USE_VTK_INSTEAD_OF_MESH
  ! VTK
  character(len=MAX_STRING_LEN) :: mesh_file
  ! global point data
  real,dimension(:),allocatable :: total_dat
  real,dimension(:,:),allocatable :: total_dat_xyz
  integer,dimension(:,:),allocatable :: total_dat_con
#else
  !!! .mesh specific !!!!!!!!!!!
  integer :: pfd, efd
  character(len=MAX_STRING_LEN) :: command_name
  character(len=MAX_STRING_LEN) :: pt_mesh_file1, pt_mesh_file2, mesh_file, em_mesh_file
#endif

#ifdef ADIOS_INPUT
  integer :: sizeprocs
  character(len=MAX_STRING_LEN) :: value_file_name, mesh_file_name
  integer(kind=8) :: value_handle, mesh_handle
#else
  integer :: iregion,njunk
  character(len=MAX_STRING_LEN) :: prname_topo, prname_file, dimension_file
  character(len=MAX_STRING_LEN) :: in_file_dir, in_topo_dir
  character(len=MAX_STRING_LEN) :: sline,slice_list_name
#endif

  ! starts here---------------------------------------------------------------
  ier = 0 ! avoids compiler warning in case of ADIOS and VTK output

  ! ADIOS MPI initialization
#ifdef ADIOS_INPUT
  ! starts mpi
  call init_mpi()
  call world_size(sizeprocs)
  ! checks number of processes
  ! note: must run as a single process with: mpirun -np 1 ..
  if (sizeprocs /= 1) then
    call world_rank(myrank)
    ! usage info
    if (myrank == 0) then
      print *, "ADIOS requires MPI functionality. However, this program executes as sequential program."
      print *, "Invalid number of processes used: ", sizeprocs, " procs"
      print *
      print *, "Please run: mpirun -np 1 ./bin/xcombine_vol_data_**_adios"
    endif
    call abort_mpi()
  endif
#endif

  ! checks array sizes
  if (NSPEC_CRUST_MANTLE < NSPEC_OUTER_CORE .or. NSPEC_CRUST_MANTLE < NSPEC_INNER_CORE) &
             stop 'This program needs that NSPEC_CRUST_MANTLE > NSPEC_OUTER_CORE and NSPEC_INNER_CORE'

  ! reads input arguments
#ifndef ADIOS_INPUT
  do i = 1, 7
    call get_command_argument(i,arg(i))

    ! usage info
    if (i < 7 .and. len_trim(arg(i)) == 0) then
      print *, ' '
      print *, ' Usage: xcombine_vol_data slice_list filename input_topo_dir input_file_dir '
      print *, '        output_dir high/low-resolution [region]'
      print *, ' with'
      print *, '   slice_list          - text file containing slice numbers to combine (or use name "all" for all slices)'
      print *, '   filename            - root file name of files proc***_filename.bin ("vp", "vsh", "alpha_kernel",..)'
      print *, '   input_topo_dir      - directory containing mesh topology (e.g. DATABASES_MPI/)'
      print *, '   input_file_dir      - directory containing data files to combine (e.g. also in directory DATABASES_MPI/)'
      print *, '   output_dir          - directory for output files (e.g. OUTPUT_FILES/)'
      print *, '   high/low-resolution - 0 == low, 1 == high-resolution'
      print *, '   region              - (optional) region number, only use 1 == crust/mantle, 2 == outer core, 3 == inner core'
      print *, ' '
      print *, ' ***** Notice: now allow different input dir for topo and kernel files ******** '
      print *, '   expect to have the topology and filename.bin(NGLLX,NGLLY,NGLLZ,nspec) '
      print *, '   already collected to input_topo_dir and input_file_dir'
      print *, '   output mesh files (filename_points.mesh, filename_elements.mesh) are saved to output_dir '
      print *, '   give 0 for low resolution and 1 for high resolution'
      print *, '   if region is not specified, all 3 regions will be collected, otherwise, only collect regions specified'
      stop ' Reenter command line options'
    endif
  enddo

  ! get region id
  if (len_trim(arg(7)) == 0) then
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
  slice_list_name = trim(arg(1))
  if (trim(slice_list_name) == 'all') then
    ! uses all slices to combine
    do i = 0,NPROCTOT_VAL-1
      num_node = num_node + 1
      if (num_node > MAX_NUM_NODES ) stop 'Error number of slices exceeds MAX_NUM_NODES...'
      node_list(num_node) = i
    enddo
  else
    ! reads in slice file with specified slices
    open(unit = IIN, file = trim(slice_list_name), status = 'old',iostat = ier)
    if (ier /= 0) then
      print *,'no file: ',trim(slice_list_name)
      stop 'Error opening slices file'
    endif

    ! file format uses one slice number per line
    ! example content of 'slices.txt':
    !  0
    !  11
    !  14
    ! end of file to combine 3 slices (from process numbers 0, 11 and 14)
    do while (1 == 1)
      read(IIN,'(a)',iostat=ier) sline
      if (ier /= 0) exit
      read(sline,*,iostat=ier) njunk
      if (ier /= 0) exit
      num_node = num_node + 1
      if (num_node > MAX_NUM_NODES ) stop 'Error number of slices exceeds MAX_NUM_NODES...'
      node_list(num_node) = njunk
    enddo
    close(IIN)
  endif

  ! output info
  print *, 'slice list: '
  print *, node_list(1:num_node)
  print *, ' '

  ! file to collect
  var_name = arg(2)

  ! input and output dir
  in_topo_dir= arg(3)
  in_file_dir= arg(4)
  outdir = arg(5)

  ! resolution
  read(arg(6),*) ires
#else
  ! ADIOS input arguments
  do i = 1, 7
    call get_command_argument(i,arg(i))
  enddo
  call read_args_adios(arg, MAX_NUM_NODES, node_list, num_node, &
                       var_name, value_file_name, mesh_file_name, &
                       outdir, ires, irs, ire, NPROCTOT_VAL)
#endif
  filename = var_name

  ! output info
  print *, 'regions: start =', irs, ' to end =', ire

  ! sets up loop increments
  di = 0
  dj = 0
  dk = 0
  if (ires == 0) then
    print *, 'using mesh with: low resolution'
    HIGH_RESOLUTION_MESH = .false.
    di = NGLLX-1
    dj = NGLLY-1
    dk = NGLLZ-1
  else if (ires == 1) then
    print *, 'using mesh with: high resolution'
    HIGH_RESOLUTION_MESH = .true.
    di = 1
    dj = 1
    dk = 1
  else if (ires == 2) then
    print *, 'using mesh with: mid resolution'
    HIGH_RESOLUTION_MESH = .false.
    di = int((NGLLX-1)/2.0)
    dj = int((NGLLY-1)/2.0)
    dk = int((NGLLZ-1)/2.0)
  else
    stop 'resolution setting must be 0, 1, or 2'
  endif

  ! sets up ellipticity splines in order to remove ellipticity from point coordinates
  if (CORRECT_ELLIPTICITY) call make_ellipticity(nspl,rspl,espl,espl2,ONE_CRUST)

#ifdef ADIOS_INPUT
  call init_adios(value_file_name, mesh_file_name, value_handle, mesh_handle)
#endif

  do ir = irs, ire
    print *, '----------- Region ', ir, '----------------'

#ifdef USE_VTK_INSTEAD_OF_MESH
    ! VTK
    ! not special file names required
    !continue
#else
    !!! .mesh specific !!!!!!!!!!!
    ! open Paraview output mesh file
    write(pt_mesh_file1,'(a,i1,a)') trim(outdir)//'/' // 'reg_',ir,'_'//trim(filename)//'_point1.mesh'
    write(pt_mesh_file2,'(a,i1,a)') trim(outdir)//'/' // 'reg_',ir,'_'//trim(filename)//'_point2.mesh'
    write(mesh_file,'(a,i1,a)') trim(outdir)//'/' // 'reg_',ir,'_'//trim(filename)//'.mesh'
    write(em_mesh_file,'(a,i1,a)') trim(outdir)//'/' // 'reg_',ir,'_'//trim(filename)//'_element.mesh'

    call open_file_fd(trim(pt_mesh_file1)//char(0),pfd)
    call open_file_fd(trim(em_mesh_file)//char(0),efd)
#endif

    ! figure out total number of points and elements for high-res mesh

    do it = 1, num_node

      iproc = node_list(it)

      print *, 'Reading slice ', iproc
#ifndef ADIOS_INPUT
      write(prname_topo,'(a,i6.6,a,i1,a)') trim(in_topo_dir)//'/proc',iproc,'_reg',ir,'_'
      write(prname_file,'(a,i6.6,a,i1,a)') trim(in_file_dir)//'/proc',iproc,'_reg',ir,'_'

      dimension_file = trim(prname_topo) //'solver_data.bin'
      open(unit = IIN,file = trim(dimension_file),status='old',action='read', iostat = ier, form='unformatted')
      if (ier /= 0) then
       print *,'Error ',ier
       print *,'file:',trim(dimension_file)
       stop 'Error opening file'
      endif
      read(IIN) nspec_list(it)
      read(IIN) nglob_list(it)
      close(IIN)
#else
      ! adios reading
      call read_scalars_adios_mesh(mesh_handle, iproc, ir, nglob_list(it), nspec_list(it))
#endif

      ! check
      if (nspec_list(it) > NSPEC_CRUST_MANTLE ) stop 'Error file nspec too big, please check compilation'
      if (nglob_list(it) > NGLOB_CRUST_MANTLE ) stop 'Error file nglob too big, please check compilation'

      ! theoretical number of points and elements
      npoint(it) = nglob_list(it)
      if (ires == 0) then
        ! low resolution
        nelement(it) = nspec_list(it)
      else if (ires == 1) then
        ! high resolution
        nelement(it) = nspec_list(it) * (NGLLX-1) * (NGLLY-1) * (NGLLZ-1)
      else if (ires == 2) then
        ! mid resolution
        nelement(it) = nspec_list(it) * (NGLLX-1) * (NGLLY-1) * (NGLLZ-1) / 8
      endif

    enddo

    print *, 'nspec_list(it) = ', nspec_list(1:num_node)
    print *, 'nglob_list(it) = ', nglob_list(1:num_node)
    print *, 'nelement(it)   = ', nelement(1:num_node)

#ifdef USE_VTK_INSTEAD_OF_MESH
    ! VTK
    print *
    print *,'VTK initial total points: ',sum(npoint(1:num_node))
    print *,'VTK initial total elements: ',sum(nelement(1:num_node))
    print *
    ! creates array to hold point data
    allocate(total_dat(sum(npoint(1:num_node))),stat=ier)
    if (ier /= 0 ) stop 'Error allocating total_dat array'
    total_dat(:) = 0.0
    allocate(total_dat_xyz(3,sum(npoint(1:num_node))),stat=ier)
    if (ier /= 0 ) stop 'Error allocating total_dat_xyz array'
    total_dat_xyz(:,:) = 0.0
    allocate(total_dat_con(8,sum(nelement(1:num_node))),stat=ier)
    if (ier /= 0 ) stop 'Error allocating total_dat_con array'
    total_dat_con(:,:) = 0
#else
    ! .mesh
    call write_integer_fd(efd,sum(nelement(1:num_node)))
#endif

    np = 0
    ne = 0

    ! write points information
    do it = 1, num_node
      ! gets slice id (process id)
      iproc = node_list(it)

      ! initializes data values
      data(:,:,:,:) = -1.e10

      ! output info
      print *, ' '
      print *, 'Reading slice ', iproc

      ! reads in kernel/data values
#ifndef ADIOS_INPUT
      write(prname_topo,'(a,i6.6,a,i1,a)') trim(in_topo_dir)//'/proc',iproc,'_reg',ir,'_'
      write(prname_file,'(a,i6.6,a,i1,a)') trim(in_file_dir)//'/proc',iproc,'_reg',ir,'_'

      ! filename.bin
      data_file = trim(prname_file) // trim(filename) // '.bin'

      open(unit = IIN,file = trim(data_file),status='old',action='read', iostat = ier,form ='unformatted')
      if (ier /= 0) then
        print *,'Error opening file: ',trim(data_file)
        stop 'Error opening file'
      endif
      read(IIN,iostat=ier) data(:,:,:,1:nspec_list(it))
      if (ier /= 0) then
        print *,'Error reading file: ',trim(data_file)
        stop 'Error reading data'
      endif
      close(IIN)
#else
      ! adios reading
      call read_values_adios(value_handle, var_name, iproc, ir, nspec_list(it), data)
      data_file = trim(var_name)
#endif
      ! output info
      print *,trim(data_file)
      print *,'  min/max value: ',minval(data(:,:,:,1:nspec_list(it))),maxval(data(:,:,:,1:nspec_list(it)))
      print *

      ! topology file
      ! reads in mesh coordinates and local-to-global mapping (ibool)
#ifndef ADIOS_INPUT
      data_file = trim(prname_topo) // 'solver_data.bin'
      !print *, trim(data_file)

      open(unit = IIN,file = trim(data_file),status='old',action='read', iostat = ier, form='unformatted')
      if (ier /= 0) then
        print *,'Error opening file: ',trim(data_file)
        stop 'Error opening topo file'
      endif
      xstore(:) = 0.0
      ystore(:) = 0.0
      zstore(:) = 0.0
      ibool(:,:,:,:) = -1
      read(IIN) nspec_list(it)
      read(IIN) nglob_list(it)
      read(IIN) xstore(1:nglob_list(it))
      read(IIN) ystore(1:nglob_list(it))
      read(IIN) zstore(1:nglob_list(it))
      read(IIN) ibool(:,:,:,1:nspec_list(it))
      if (ir == 3) read(IIN) idoubling_inner_core(1:nspec_list(it)) ! flag that can indicate fictitious elements
      close(IIN)
#else
      ! adios reading
      call read_coordinates_adios_mesh(mesh_handle, iproc, ir, nglob_list(it), nspec_list(it), xstore, ystore, zstore, ibool)
#endif

      !average data on global points
      if (AVERAGE_GLOBALPOINTS) then
        print *,'  averaging data on global points...'

        ! initializes averaged data
        if (.not. allocated(data_avg)) then
          allocate(data_avg(NGLOB_CRUST_MANTLE), &
                   ibool_count(NGLOB_CRUST_MANTLE),stat=ier)
          if (ier /= 0) stop 'Error allocating data_avg, ibool_count arrays'
        endif
        data_avg(:) = 0.0
        ibool_count(:) = 0

        do ispec = 1,nspec_list(it)
          ! checks if element counts
          if (ir == 3) then
            ! inner core
            ! nothing to do for fictitious elements in central cube
            if (idoubling_inner_core(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle
          endif
          ! counts and sums global point data
          do k = 1, NGLLZ, dk
            do j = 1, NGLLY, dj
              do i = 1, NGLLX, di
                iglob = ibool(i,j,k,ispec)

                dat = data(i,j,k,ispec)

                data_avg(iglob) = data_avg(iglob) + dat
                ibool_count(iglob) = ibool_count(iglob) + 1
              enddo
            enddo
          enddo
        enddo
        ! averages data
        do iglob = 1,nglob_List(it)
          if (ibool_count(iglob) > 0) then
            data_avg(iglob) = data_avg(iglob)/ibool_count(iglob)
          endif
        enddo
      endif

      mask_ibool(:) = .false.
      num_ibool(:) = 0
      numpoin = 0

      ! write point file
      do ispec = 1,nspec_list(it)
        ! checks if element counts
        if (ir == 3) then
          ! inner core
          ! nothing to do for fictitious elements in central cube
          if (idoubling_inner_core(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle
        endif

        ! writes out global point data
        do k = 1, NGLLZ, dk
          do j = 1, NGLLY, dj
            do i = 1, NGLLX, di
              iglob = ibool(i,j,k,ispec)
              if (iglob == -1 ) cycle

              if (.not. mask_ibool(iglob)) then
                numpoin = numpoin + 1
                x = xstore(iglob)
                y = ystore(iglob)
                z = zstore(iglob)

                ! remove ellipticity
                if (CORRECT_ELLIPTICITY) call revert_ellipticity(x,y,z,nspl,rspl,espl,espl2)

                ! takes the averaged data value for mesh
                if (AVERAGE_GLOBALPOINTS) then
                  dat = data_avg(iglob)
                else
                  dat = data(i,j,k,ispec)
                endif

#ifdef USE_VTK_INSTEAD_OF_MESH
                ! VTK
                total_dat(np+numpoin) = dat
                total_dat_xyz(1,np+numpoin) = x
                total_dat_xyz(2,np+numpoin) = y
                total_dat_xyz(3,np+numpoin) = z
#else
                ! .mesh
                call write_real_fd(pfd,x)
                call write_real_fd(pfd,y)
                call write_real_fd(pfd,z)
                call write_real_fd(pfd,dat)
#endif
                mask_ibool(iglob) = .true.
                num_ibool(iglob) = numpoin
              endif

            enddo ! i
          enddo ! j
        enddo ! k
      enddo !ispec

      ! no way to check the number of points for low-res
      if (HIGH_RESOLUTION_MESH) then
        if (ir == 3) then
          npoint(it) = numpoin
        else if (numpoin /= npoint(it)) then
          print *,'region:',ir
          print *,'Error number of points:',numpoin,npoint(it)
          stop 'Error different number of points (high-res)'
        endif
      else
        npoint(it) = numpoin
      endif

      ! write elements file
      numpoin = 0
      do ispec = 1, nspec_list(it)
        ! checks if element counts
        if (ir == 3) then
          ! inner core
          ! fictitious elements in central cube
          if (idoubling_inner_core(ispec) == IFLAG_IN_FICTITIOUS_CUBE) then
            ! connectivity must be given, otherwise element count would be wrong
            ! maps "fictitious" connectivity, element is all with iglob = 1

#ifdef USE_VTK_INSTEAD_OF_MESH
            ! VTK, can avoid output of ficticious inner core elements
#else
            ! .mesh
            do k = 1, NGLLZ-1, dk
              do j = 1, NGLLY-1, dj
                do i = 1, NGLLX-1, di
                  call write_integer_fd(efd,1)
                  call write_integer_fd(efd,1)
                  call write_integer_fd(efd,1)
                  call write_integer_fd(efd,1)
                  call write_integer_fd(efd,1)
                  call write_integer_fd(efd,1)
                  call write_integer_fd(efd,1)
                  call write_integer_fd(efd,1)
                enddo ! i
              enddo ! j
            enddo ! k
#endif
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

#ifdef USE_VTK_INSTEAD_OF_MESH
              ! VTK
              ! note: indices for VTK start at 0
              total_dat_con(1,numpoin + ne) = n1
              total_dat_con(2,numpoin + ne) = n2
              total_dat_con(3,numpoin + ne) = n3
              total_dat_con(4,numpoin + ne) = n4
              total_dat_con(5,numpoin + ne) = n5
              total_dat_con(6,numpoin + ne) = n6
              total_dat_con(7,numpoin + ne) = n7
              total_dat_con(8,numpoin + ne) = n8
#else
              ! .mesh
              call write_integer_fd(efd,n1)
              call write_integer_fd(efd,n2)
              call write_integer_fd(efd,n3)
              call write_integer_fd(efd,n4)
              call write_integer_fd(efd,n5)
              call write_integer_fd(efd,n6)
              call write_integer_fd(efd,n7)
              call write_integer_fd(efd,n8)
#endif
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

#ifdef USE_VTK_INSTEAD_OF_MESH
    ! VTK
    ! outputs unstructured grid file
    write(mesh_file,'(a,i1,a)') trim(outdir)//'/' // 'reg_',ir,'_'//trim(filename)//'.vtk'

    call write_VTK_movie_data(ne,np,total_dat_xyz,total_dat_con,total_dat,mesh_file,var_name)

    ! free arrays for this region
    deallocate(total_dat,total_dat_xyz,total_dat_con)

    print *,'written: ',trim(mesh_file)
    print *
#else
    call close_file_fd(pfd)
    call close_file_fd(efd)

    ! add the critical piece: total number of points
    call open_file_fd(trim(pt_mesh_file2)//char(0),pfd)
    call write_integer_fd(pfd,np)
    call close_file_fd(pfd)

    command_name='cat '//trim(pt_mesh_file2)//' '//trim(pt_mesh_file1)//' '//trim(em_mesh_file)//'>'//trim(mesh_file)
    print *, ' '
    print *, 'cat mesh files: '
    print *, trim(command_name)
    call system_command(command_name)
#endif
  enddo

#ifdef ADIOS_INPUT
  call clean_adios(value_handle, mesh_handle)
  ! shuts down mpi
  call finalize_mpi()
#endif

  print *, 'Done writing mesh files'
  print *, ' '

end program combine_vol_data
