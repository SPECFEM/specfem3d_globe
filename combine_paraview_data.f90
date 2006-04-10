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

program combine_paraview_data

! combines the database files on several slices.
! the local database file needs to have been collected onto the frontend (copy_local_database.pl)

  implicit none

  include 'constants.h'
  include 'OUTPUT_FILES/values_from_mesher.h'

  integer i,j,k,ispec, ios, it
  integer iproc, num_node, node_list(300), nspec(300), nglob(300), nglob1(300),npoint_all, nelement_all
  integer np, ne, npoint(300), nelement(300), njunk, njunk2, n1, n2, n3, n4, n5, n6, n7, n8
  integer ibool(NGLLX,NGLLY,NGLLZ,NSPECMAX_CRUST_MANTLE)

  integer numpoin, iglob1, iglob2, iglob3, iglob4, iglob5, iglob6, iglob7, iglob8, iglob
  logical mask_ibool(NGLOBMAX_CRUST_MANTLE)
  real(kind=CUSTOM_REAL) data(NGLLX,NGLLY,NGLLZ,NSPECMAX_CRUST_MANTLE)
  real(kind=CUSTOM_REAL),dimension(NGLOBMAX_CRUST_MANTLE) :: xstore, ystore, zstore
  real x, y, z, dat(NGLLX,NGLLY,NGLLZ,NSPECMAX_CRUST_MANTLE)
  character(len=150) :: sline, arg(6), filename, indir, outdir, prname, dimension_file
  character(len=150) :: mesh_file, local_point_file, local_element_file, local_file, local_data_file, local_ibool_file
  integer :: num_ibool(NGLOBMAX_CRUST_MANTLE)
  logical :: HIGH_RESOLUTION_MESH
  integer :: ires, iregion,irs,ire,ir


  do i = 1, 6
    call getarg(i,arg(i))
    if (i < 6 .and. trim(arg(i)) == '') then
      print *, ' '
      print *, ' Usage: xcombine_data slice_list filename input_dir output_dir high/low-resolution [region]'
      print *, ' possible filenames -- rho_vp, rho_vs, kappastore, mustore etc'
      print *, '   stored in the local directory as real(kind=CUSTOM_REAL) filename(NGLLX,NGLLY,NGLLZ,nspec)  '
      print *, '   in filename.bin'
      print *, ' files have been collected in input_dir, output mesh file goes to output_dir '
      print *, ' give 0 for low resolution and 1 for high resolution'
      print *, ' if region is not specified, all 3 regions will be collected, otherwise, only collect regions specified'
      stop ' Reenter command line options'
    endif
  enddo

  if (NSPECMAX_CRUST_MANTLE < NSPECMAX_OUTER_CORE .or. NSPECMAX_CRUST_MANTLE < NSPEC_INNER_CORE) &
    stop 'This program needs that NSPECMAX_CRUST_MANTLE > NSPECMAX_OUTER_CORE and NSPEC_INNER_CORE'

! get slice list
  if (trim(arg(6)) == '') then
    iregion  = 0
  else
    read(arg(6),*) iregion
  endif
  if (iregion > 3 .or. iregion < 0) stop 'Iregion = 0,1,2,3'
  num_node = 0
  open(unit = 20, file = trim(arg(1)), status = 'unknown',iostat = ios)
  if (ios /= 0) stop 'Error opening file'
  do while ( 1 == 1)
    read(20,'(a)',iostat=ios) sline
    if (ios /= 0) exit
    read(sline,*,iostat=ios) njunk
    if (ios /= 0) exit
    num_node = num_node + 1
    node_list(num_node) = njunk
  enddo
  close(20)
  filename = arg(2)
  indir= arg(3)
  outdir = arg(4)
  read(arg(5),*) ires

  if (ires == 0) then
    HIGH_RESOLUTION_MESH = .false.
  else
    HIGH_RESOLUTION_MESH = .true.
  endif

  print *, 'Slice list: '
  print *, node_list(1:num_node)
  print *, ' '

  ! write point and scalar information

  if (iregion == 0) then
    irs = 1; ire = 3
  else
    irs = iregion; ire = irs
  endif

  do ir = irs, ire
    print *, '----------- Region ', ir, '----------------'

  ! open paraview output mesh file
    write(mesh_file,'(a,i1,a)') trim(outdir)//'/' // 'reg_',ir,'_'//trim(filename)//'.mesh'
    call open_file(trim(mesh_file)//char(0))


  ! figure out total number of points

  do it = 1, num_node

    iproc = node_list(it)

    print *, 'Reading slice ', iproc
    write(prname,'(a,i4.4,a,i1,a)') trim(indir)//'/proc',iproc,'_reg',ir,'_'

    dimension_file = trim(prname) //'array_dims.txt'
    open(unit = 27,file = trim(dimension_file),status='old', iostat = ios)
    if (ios /= 0) stop 'Error opening file'

    read(27,*) nspec(it)
    read(27,*) nglob(it)
    if (ir == 3) then
      read(27,*) nglob1(it)
      if (.not. HIGH_RESOLUTION_MESH) nspec(it) = nspec_inner_core
    else
      nglob1(it) = nglob(it)
    endif
    close(27)

    if (.not. HIGH_RESOLUTION_MESH) then
      local_point_file = trim(prname) //'AVS_DXpoints' // '.txt'
      open(unit = 27,file = trim(local_point_file),status='old', iostat = ios)
      if (ios /= 0) stop 'Error opening file'

      read(27,*) npoint(it)
      close(27)
    else
      npoint(it) = nglob1(it)
    endif

  enddo

  npoint_all = sum(npoint(1:num_node))
  print *, 'Region ', ir, ' total number of points = ', npoint_all
  print *, 'nspec(it) = ', nspec(1:num_node)
  print *, 'nglob(it) = ', nglob1(1:num_node)

  np = 0

  ! write point and scalar information

  do it = 1, num_node

    iproc = node_list(it)

    print *, ' '
    print *, 'Reading slice ', iproc
    write(prname,'(a,i4.4,a,i1,a)') trim(indir)//'/proc',iproc,'_reg',ir,'_'

    local_data_file = trim(prname) // trim(filename) // '.bin'
    open(unit = 27,file = trim(local_data_file),status='old', iostat = ios,form ='unformatted')
    if (ios /= 0) stop 'Error opening file'

    read(27) data(:,:,:,1:nspec(it))
    close(27)
    print *, trim(local_data_file)

    dat = data

  ! ibool file
    local_ibool_file = trim(prname) // 'ibool' // '.bin'
    open(unit = 28,file = trim(local_ibool_file),status='old', iostat = ios, form='unformatted')
    if (ios /= 0) stop 'Error opening file'

    read(28) ibool(:,:,:,1:nspec(it))
    close(28)
    print *, trim(local_ibool_file)

    mask_ibool(:) = .false.
    numpoin = 0

    if (.not. HIGH_RESOLUTION_MESH) then

      local_point_file = trim(prname) // 'AVS_DXpoints.txt'
      open(unit = 25, file = trim(local_point_file), status = 'old', iostat = ios)
      if (ios /= 0) stop 'Error opening file'

      read(25,*) njunk

      if (it == 1) then
        call write_integer(npoint_all)
      endif

      do ispec=1,nspec(it)
        iglob1=ibool(1,1,1,ispec)
        iglob2=ibool(NGLLX,1,1,ispec)
        iglob3=ibool(NGLLX,NGLLY,1,ispec)
        iglob4=ibool(1,NGLLY,1,ispec)
        iglob5=ibool(1,1,NGLLZ,ispec)
        iglob6=ibool(NGLLX,1,NGLLZ,ispec)
        iglob7=ibool(NGLLX,NGLLY,NGLLZ,ispec)
        iglob8=ibool(1,NGLLY,NGLLZ,ispec)

        if(.not. mask_ibool(iglob1)) then
          numpoin = numpoin + 1
          read(25,*) njunk, x, y, z
          call write_real(x)
          call write_real(y)
          call write_real(z)
          call write_real(dat(1,1,1,ispec))
          mask_ibool(iglob1) = .true.
        endif
        if(.not. mask_ibool(iglob2)) then
          numpoin = numpoin + 1
          read(25,*) njunk, x, y, z
          call write_real(x)
          call write_real(y)
          call write_real(z)
          call write_real(dat(NGLLX,1,1,ispec))
          mask_ibool(iglob2) = .true.
        endif
        if(.not. mask_ibool(iglob3)) then
          numpoin = numpoin + 1
          read(25,*) njunk, x, y, z
          call write_real(x)
          call write_real(y)
          call write_real(z)
          call write_real(dat(NGLLX,NGLLY,1,ispec))
          mask_ibool(iglob3) = .true.
        endif
        if(.not. mask_ibool(iglob4)) then
          numpoin = numpoin + 1
          read(25,*) njunk, x, y, z
          call write_real(x)
          call write_real(y)
          call write_real(z)
          call write_real(dat(1,NGLLY,1,ispec))
          mask_ibool(iglob4) = .true.
        endif
        if(.not. mask_ibool(iglob5)) then
          numpoin = numpoin + 1
          read(25,*) njunk, x, y, z
          call write_real(x)
          call write_real(y)
          call write_real(z)
          call write_real(dat(1,1,NGLLZ,ispec))
          mask_ibool(iglob5) = .true.
        endif
        if(.not. mask_ibool(iglob6)) then
          numpoin = numpoin + 1
          read(25,*) njunk, x, y, z
          call write_real(x)
          call write_real(y)
          call write_real(z)
          call write_real(dat(NGLLX,1,NGLLZ,ispec))
          mask_ibool(iglob6) = .true.
        endif
        if(.not. mask_ibool(iglob7)) then
          numpoin = numpoin + 1
          read(25,*) njunk, x, y, z
          call write_real(x)
          call write_real(y)
          call write_real(z)
          call write_real(dat(NGLLX,NGLLY,NGLLZ,ispec))
          mask_ibool(iglob7) = .true.
        endif
        if(.not. mask_ibool(iglob8)) then
          numpoin = numpoin + 1
          read(25,*) njunk, x, y, z
          call write_real(x)
          call write_real(y)
          call write_real(z)
          call write_real(dat(1,NGLLY,NGLLZ,ispec))
          mask_ibool(iglob8) = .true.
        endif
      enddo ! ispec
      close(25)

    else  ! high resolution

      if (it == 1) then
        call write_integer(npoint_all)
      endif

      local_file = trim(prname)//'x.bin'
      open(unit = 27,file = trim(prname)//'x.bin',status='old', iostat = ios,form ='unformatted')
     if (ios /= 0) stop 'Error opening file'

      read(27) xstore(1:nglob(it))
      close(27)
      local_file = trim(prname)//'y.bin'
      open(unit = 27,file = trim(prname)//'y.bin',status='old', iostat = ios,form ='unformatted')
     if (ios /= 0) stop 'Error opening file'

      read(27) ystore(1:nglob(it))
      close(27)
      local_file = trim(prname)//'z.bin'
      open(unit = 27,file = trim(prname)//'z.bin',status='old', iostat = ios,form ='unformatted')
      if (ios /= 0) stop 'Error opening file'

      read(27) zstore(1:nglob(it))
      close(27)

      do ispec=1,nspec(it)
        do k = 1, NGLLZ
          do j = 1, NGLLY
            do i = 1, NGLLX
              iglob = ibool(i,j,k,ispec)
              if(.not. mask_ibool(iglob)) then
                numpoin = numpoin + 1
                x = xstore(iglob)
                y = ystore(iglob)
                z = zstore(iglob)
                call write_real(x)
                call write_real(y)
                call write_real(z)
                call write_real(dat(i,j,k,ispec))
                mask_ibool(iglob) = .true.
              endif
            enddo ! i
          enddo ! j
        enddo ! k
      enddo !ispec
    endif

    if (numpoin /= npoint(it)) stop 'different number of points'
    np = np + npoint(it)

  enddo  ! all slices for points

 if (np /=  npoint_all) stop 'Error: Number of total points are not consistent'
 print *, 'Total number of points: ', np
 print *, ' '

! figure out total number of elements

  do it = 1, num_node

    iproc = node_list(it)

    print *, 'Reading slice ', iproc
    write(prname,'(a,i4.4,a,i1,a)') trim(indir)//'/proc',iproc,'_reg',ir,'_'

    if (.not. HIGH_RESOLUTION_MESH) then
      local_point_file = trim(prname) //'AVS_DXelements' // '.txt'
      open(unit = 27,file = trim(local_point_file),status='old', iostat = ios)
     if (ios /= 0) stop 'Error opening file'

      read(27,*) nelement(it)
      nelement(it) = nspec(it)
      close(27)
    else
      nelement(it) = nspec(it) * (NGLLX-1) * (NGLLY-1) * (NGLLZ-1)
    endif

  enddo
  nelement_all = sum(nelement(1:num_node))

  ne = 0
! write element information
 do it = 1, num_node

    iproc = node_list(it)

    print *, 'Reading slice ', iproc
    write(prname,'(a,i4.4,a,i1,a)') trim(indir)//'/proc',iproc,'_reg',ir,'_'

    if (it == 1) then
      np = 0
    else
      np = sum(npoint(1:it-1))
    endif

    if (.not. HIGH_RESOLUTION_MESH) then

      local_element_file = trim(prname) // 'AVS_DXelements.txt'
      open(unit = 26, file = trim(local_element_file), status = 'old', iostat = ios)
      if (ios /= 0) stop 'Error opening file'

      print *, trim(local_element_file)

      read(26, *) njunk
      if (it == 1) then
        call write_integer(nelement_all)
      end if

      do i = 1, nelement(it)
        read(26,*) njunk, njunk2, n1, n2, n3, n4, n5, n6, n7, n8
        n1 = n1+np-1; n2 = n2+np-1; n3 = n3+np-1; n4 = n4+np-1
        n5 = n5+np-1; n6 = n6+np-1; n7 = n7+np-1; n8 = n8+np-1
        call write_integer(n1)
        call write_integer(n2)
        call write_integer(n3)
        call write_integer(n4)
        call write_integer(n5)
        call write_integer(n6)
        call write_integer(n7)
        call write_integer(n8)
      enddo
      close(26)

    else ! high resolution mesh

      if (it == 1) then
        call write_integer(nelement_all)
      endif

      local_ibool_file = trim(prname) // 'ibool' // '.bin'
      open(unit = 28,file = trim(local_ibool_file),status='old', iostat = ios, form='unformatted')
       if (ios /= 0) stop 'Error opening file'

      read(28) ibool(:,:,:,1:nspec(it))
      close(28)

      numpoin = 0
      mask_ibool = .false.
      do ispec=1,nspec(it)
        do k = 1, NGLLZ
          do j = 1, NGLLY
            do i = 1, NGLLX
              iglob = ibool(i,j,k,ispec)
              if(.not. mask_ibool(iglob)) then
                numpoin = numpoin + 1
                num_ibool(iglob) = numpoin
                mask_ibool(iglob) = .true.
              endif
            enddo ! i
          enddo ! j
        enddo ! k
      enddo !ispec

      do ispec = 1, nspec(it)
        do k = 1, NGLLZ-1
          do j = 1, NGLLY-1
            do i = 1, NGLLX-1
              iglob1 = ibool(i,j,k,ispec)
              iglob2 = ibool(i+1,j,k,ispec)
              iglob3 = ibool(i+1,j+1,k,ispec)
              iglob4 = ibool(i,j+1,k,ispec)
              iglob5 = ibool(i,j,k+1,ispec)
              iglob6 = ibool(i+1,j,k+1,ispec)
              iglob7 = ibool(i+1,j+1,k+1,ispec)
              iglob8 = ibool(i,j+1,k+1,ispec)
              n1 = num_ibool(iglob1)+np-1
              n2 = num_ibool(iglob2)+np-1
              n3 = num_ibool(iglob3)+np-1
              n4 = num_ibool(iglob4)+np-1
              n5 = num_ibool(iglob5)+np-1
              n6 = num_ibool(iglob6)+np-1
              n7 = num_ibool(iglob7)+np-1
              n8 = num_ibool(iglob8)+np-1
              call write_integer(n1)
              call write_integer(n2)
              call write_integer(n3)
              call write_integer(n4)
              call write_integer(n5)
              call write_integer(n6)
              call write_integer(n7)
              call write_integer(n8)
            enddo
          enddo
        enddo
      enddo

    endif
    ne = ne + nelement(it)

  enddo ! num_node
  if (ne /= nelement_all) stop 'Number of total elements are not consistent'
  print *, 'Total number of elements: ', ne

  call close_file()

  print *, 'Done writing '//trim(mesh_file)
  print *, ' '

  enddo
  print *, ' '

end program combine_paraview_data

