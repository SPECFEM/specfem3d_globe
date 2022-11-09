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


  subroutine prepare_vtk_window()

! prepares arrays for VTK run-time visualization

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_movie

  implicit none

  ! local parameters
  integer :: i,j,k,iglob,ispec,inum,ier
  integer :: id1,id2,id3,id4,id5,id6,id7,id8
  integer :: ispec2D,NIT_res

  ! free surface points
  integer :: free_np,free_nspec
  real, dimension(:),allocatable :: free_x,free_y,free_z
  integer, dimension(:,:),allocatable :: free_conn
  integer, dimension(:),allocatable :: free_perm
  ! gather arrays for multi-MPI simulations
  real, dimension(:),allocatable :: free_x_all,free_y_all,free_z_all
  integer, dimension(:,:),allocatable :: free_conn_all
  integer, dimension(:),allocatable :: free_conn_offset_all,free_conn_nspec_all
  integer, dimension(:),allocatable :: free_points_all,free_offset_all
  integer :: free_np_all,free_nspec_all

  ! volume points
  integer :: vol_np,vol_nspec
  real, dimension(:),allocatable :: vol_x,vol_y,vol_z
  integer, dimension(:,:),allocatable :: vol_conn
  integer, dimension(:),allocatable :: vol_perm
  ! gather arrays for multi-MPI simulations
  real, dimension(:),allocatable :: vol_x_all,vol_y_all,vol_z_all
  integer, dimension(:,:),allocatable :: vol_conn_all
  integer, dimension(:),allocatable :: vol_conn_offset_all,vol_conn_nspec_all
  integer :: vol_nspec_all,ispec_start,ispec_end
  real,dimension(1) :: dummy
  integer,dimension(1) :: dummy_i

  real(kind=CUSTOM_REAL) :: x,y,z

  !-----------------------------------------------------------------------
  ! user parameter
  logical, parameter :: VTK_USE_HIRES         = .false.
  logical, parameter :: VTK_SHOW_FREESURFACE  = .true.
  logical, parameter :: VTK_SHOW_VOLUME       = .true.
  !-----------------------------------------------------------------------

  ! checks if anything to do
  if (.not. VTK_MODE) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "preparing VTK runtime visualization"
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! to avoid compiler warnings
  !NPROC = NPROCTOT_VAL

  ! adds source
  if (myrank == 0) then
    ! user output
    write(IMAIN,*) "  VTK source sphere:"
    call prepare_vtksource(vtkdata_source_x,vtkdata_source_y,vtkdata_source_z)
  endif
  call synchronize_all()

  ! mask
  allocate(vtkmask(NGLOB_CRUST_MANTLE),stat=ier)
  if (ier /= 0 ) stop 'Error allocating arrays'

  if (VTK_USE_HIRES) then
    NIT_res = 1
  else
    NIT_res = NGLLX - 1
  endif

  ! free surface
  if (VTK_SHOW_FREESURFACE) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "  VTK free surface:"
      write(IMAIN,*) "    free surface elements    : ",NSPEC_TOP
    endif

    ! counts global free surface points
    vtkmask(:) = .false.

    ! determines number of global points on surface
    do ispec2D = 1, NSPEC_TOP ! NSPEC2D_TOP(IREGION_CRUST_MANTLE)
      ispec = ibelm_top_crust_mantle(ispec2D)
      ! in case of global, NCHUNKS_VAL == 6 simulations, be aware that for
      ! the cubed sphere, the mapping changes for different chunks,
      ! i.e. e.g. x(1,1) and x(5,5) flip left and right sides of the elements in geographical coordinates.
      ! for future consideration, like in create_movie_GMT_global.f90 ...
      k = NGLLZ
      ! loop on all the points inside the element
      do j = 1,NGLLY,NIT_res
        do i = 1,NGLLX,NIT_res
          iglob = ibool_crust_mantle(i,j,k,ispec)
          vtkmask(iglob) = .true.
        enddo
      enddo
    enddo

    ! loads free surface into data
    free_np = count(vtkmask(:))

    ! user output
    if (myrank == 0 ) write(IMAIN,*) "    loading surface points: ",free_np

    allocate(free_x(free_np),free_y(free_np),free_z(free_np),stat=ier)
    if (ier /= 0 ) stop 'Error allocating arrays'

    ! permutation array
    allocate(free_perm(NGLOB_CRUST_MANTLE),stat=ier)
    if (ier /= 0 ) stop 'Error allocating arrays'

    free_perm(:) = 0
    inum = 0
    do iglob = 1,NGLOB_CRUST_MANTLE
      if (vtkmask(iglob) .eqv. .true.) then
        inum = inum + 1
        ! note: rstore has coordinates r/theta/phi, converts back to x/y/z
        call rthetaphi_2_xyz(x,y,z,rstore_crust_mantle(1,iglob),rstore_crust_mantle(2,iglob),rstore_crust_mantle(3,iglob))
        free_x(inum) = x
        free_y(inum) = y
        free_z(inum) = z
        ! stores permutation
        free_perm(iglob) = inum
      endif
    enddo
    if (inum /= free_np) stop 'Error free_np count in loading free surface points'

    ! hi/low resolution
    if (VTK_USE_HIRES) then
      ! point connectivity
      free_nspec = NSPEC_TOP*(NGLLX-1)*(NGLLY-1)

      allocate(free_conn(4,free_nspec),stat=ier)
      if (ier /= 0 ) stop 'Error allocating arrays'

      inum = 0
      free_conn(:,:) = -1
      do ispec2D = 1,NSPEC_TOP
        ispec = ibelm_top_crust_mantle(ispec2D)
        k = NGLLZ
        do j = 1, NGLLY-1
          do i = 1, NGLLX-1
            ! indices of corner points
            id1 = free_perm(ibool_crust_mantle(i,j,k,ispec))
            id2 = free_perm(ibool_crust_mantle(i+1,j,k,ispec))
            id3 = free_perm(ibool_crust_mantle(i+1,j+1,k,ispec))
            id4 = free_perm(ibool_crust_mantle(i,j+1,k,ispec))
            ! note: indices for VTK start at 0
            inum = inum+1
            free_conn(1,inum) = id1 - 1
            free_conn(2,inum) = id2 - 1
            free_conn(3,inum) = id3 - 1
            free_conn(4,inum) = id4 - 1
          enddo
        enddo
      enddo
    else
      ! point connectivity
      free_nspec = NSPEC_TOP

      allocate(free_conn(4,free_nspec),stat=ier)
      if (ier /= 0 ) stop 'Error allocating arrays'

      inum = 0
      free_conn(:,:) = -1
      do ispec2D = 1,NSPEC_TOP
        ispec = ibelm_top_crust_mantle(ispec2D)
        ! indices of corner points
        id1 = free_perm(ibool_crust_mantle(1,1,NGLLZ,ispec))
        id2 = free_perm(ibool_crust_mantle(NGLLX,1,NGLLZ,ispec))
        id3 = free_perm(ibool_crust_mantle(NGLLX,NGLLY,NGLLZ,ispec))
        id4 = free_perm(ibool_crust_mantle(1,NGLLY,NGLLZ,ispec))
        ! note: indices for VTK start at 0
        inum = inum + 1
        free_conn(1,inum) = id1 - 1
        free_conn(2,inum) = id2 - 1
        free_conn(3,inum) = id3 - 1
        free_conn(4,inum) = id4 - 1
      enddo
    endif
    if (minval(free_conn(:,:)) < 0) stop 'Error VTK free surface point connectivity'

    ! gathers data from all MPI processes
    if (NPROC > 1) then
      ! multiple MPI processes

      ! user output
      !if (myrank == 0 ) print *,"    gathering all MPI info... "

      ! number of volume points for all partitions together
      call sum_all_i(free_np,free_np_all)
      if (myrank == 0 ) write(IMAIN,*) "    all freesurface points: ",free_np_all

      ! gathers point info
      allocate(free_points_all(NPROC),stat=ier)
      if (ier /= 0 ) stop 'Error allocating arrays'

      free_points_all(:) = 0
      call gather_all_singlei(free_np,free_points_all,NPROC)

      ! array offsets
      allocate(free_offset_all(NPROC),stat=ier)
      if (ier /= 0 ) stop 'Error allocating arrays'

      free_offset_all(1) = 0
      do i = 2, NPROC
        free_offset_all(i) = sum(free_points_all(1:i-1))
      enddo

      ! number of volume elements
      call sum_all_i(free_nspec,free_nspec_all)
      if (myrank == 0 ) write(IMAIN,*) "    all freesurface elements: ",free_nspec_all

      ! freesurface elements
      allocate(free_conn_nspec_all(NPROC),stat=ier)
      if (ier /= 0 ) stop 'Error allocating arrays'

      free_conn_nspec_all(:) = 0
      call gather_all_singlei(4*free_nspec,free_conn_nspec_all,NPROC)

      ! array offsets
      allocate(free_conn_offset_all(NPROC),stat=ier)
      if (ier /= 0 ) stop 'Error allocating arrays'

      free_conn_offset_all(1) = 0
      do i = 2, NPROC
        free_conn_offset_all(i) = sum(free_conn_nspec_all(1:i-1))
      enddo

      ! global data arrays (only needed on main process)
      if (myrank == 0) then
        ! gather locations
        allocate(free_x_all(free_np_all), &
                 free_y_all(free_np_all), &
                 free_z_all(free_np_all),stat=ier )
        if (ier /= 0 ) stop 'Error allocating free_x_all,... arrays'

        free_x_all(:) = 0.0
        free_y_all(:) = 0.0
        free_z_all(:) = 0.0

        ! connectivity
        allocate(free_conn_all(4,free_nspec_all),stat=ier)
        if (ier /= 0 ) stop 'Error allocating free_conn_all array'
        free_conn_all(:,:) = 0
      endif

      if (myrank == 0) then
        ! locations
        !if (myrank == 0 ) print *,"    locations..."
        call gatherv_all_r(free_x,free_np, &
                            free_x_all,free_points_all,free_offset_all, &
                            free_np_all,NPROC)
        call gatherv_all_r(free_y,free_np, &
                            free_y_all,free_points_all,free_offset_all, &
                            free_np_all,NPROC)
        call gatherv_all_r(free_z,free_np, &
                            free_z_all,free_points_all,free_offset_all, &
                            free_np_all,NPROC)

        ! connectivity
        !if (myrank == 0 ) print *,"    connectivity..."
        call gatherv_all_i(free_conn,4*free_nspec, &
                           free_conn_all,free_conn_nspec_all,free_conn_offset_all, &
                           free_nspec_all,NPROC)

        ! shifts connectivity ids for all additional slices
        do i = 2, NPROC
          ! divides by 4 to get nspec numbers
          ispec_start = free_conn_offset_all(i)/4 + 1
          ispec_end = free_conn_offset_all(i)/4 + free_conn_nspec_all(i)/4
          do ispec = ispec_start,ispec_end
            free_conn_all(:,ispec) = free_conn_all(:,ispec) + free_offset_all(i)
          enddo
        enddo

        !if (myrank == 0 ) print *,"    preparing VTK field..."

        ! adds free surface to VTK window
        call prepare_vtkfreesurface(free_np_all,free_x_all,free_y_all,free_z_all, &
                                    free_nspec_all,free_conn_all)

      else
        ! all other process just send data locations
        call gatherv_all_r(free_x,free_np, &
                            dummy,free_points_all,free_offset_all, &
                            1,NPROC)
        call gatherv_all_r(free_y,free_np, &
                            dummy,free_points_all,free_offset_all, &
                            1,NPROC)
        call gatherv_all_r(free_z,free_np, &
                            dummy,free_points_all,free_offset_all, &
                            1,NPROC)
        ! connectivity
        call gatherv_all_i(free_conn,4*free_nspec, &
                            dummy_i,free_conn_nspec_all,free_conn_offset_all, &
                            1,NPROC)

      endif
    else
      ! serial run
      ! creates VTK freesurface actor
      call prepare_vtkfreesurface(free_np,free_x,free_y,free_z, &
                                  free_nspec,free_conn)

    endif

    ! frees memory
    deallocate(free_x,free_y,free_z)
    deallocate(free_conn,free_perm)
    if (NPROC > 1) then
      deallocate(free_conn_nspec_all,free_conn_offset_all)
      deallocate(free_points_all,free_offset_all)
      if (myrank == 0 ) deallocate(free_x_all,free_y_all,free_z_all,free_conn_all)
    endif
  endif ! VTK_SHOW_FREESURFACE
  call synchronize_all()

  ! volume data
  if (VTK_SHOW_VOLUME) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "  VTK volume:"
      write(IMAIN,*) "    spectral elements    : ",NSPEC_CRUST_MANTLE
    endif

    ! sets new point mask
    vtkmask(:) = .false.
    do ispec = 1,NSPEC_CRUST_MANTLE
      ! hi/low resolution
      ! loops only over points
      do k = 1,NGLLZ,NIT_res
        do j = 1,NGLLY,NIT_res
          do i = 1,NGLLX,NIT_res
            iglob = ibool_crust_mantle(i,j,k,ispec)
            ! sets mask
            vtkmask(iglob) = .true.
          enddo
        enddo
      enddo
    enddo
    vol_np = count(vtkmask(:))

    ! loads volume data arrays
    if (myrank == 0 ) write(IMAIN,*) "    loading volume points: ",vol_np

    allocate(vol_x(vol_np),vol_y(vol_np),vol_z(vol_np),stat=ier)
    if (ier /= 0 ) stop 'Error allocating arrays'

    ! permutation array
    allocate(vol_perm(NGLOB_CRUST_MANTLE),stat=ier)
    if (ier /= 0 ) stop 'Error allocating arrays'

    vol_perm(:) = 0
    inum = 0
    do iglob = 1,NGLOB_CRUST_MANTLE
      if (vtkmask(iglob) .eqv. .true.) then
        inum = inum + 1
        ! note: rstore has coordinates r/theta/phi, converts back to x/y/z
        call rthetaphi_2_xyz(x,y,z,rstore_crust_mantle(1,iglob),rstore_crust_mantle(2,iglob),rstore_crust_mantle(3,iglob))
        vol_x(inum) = x
        vol_y(inum) = y
        vol_z(inum) = z
        ! stores permutation
        vol_perm(iglob) = inum
      endif
    enddo
    if (inum /= vol_np) stop 'Error vol_np count in loading volume points'

    ! hi/low resolution
    if (VTK_USE_HIRES) then
      ! point connectivity
      vol_nspec = NSPEC_CRUST_MANTLE*(NGLLX-1)*(NGLLY-1)*(NGLLZ-1)

      allocate(vol_conn(8,vol_nspec),stat=ier)
      if (ier /= 0 ) stop 'Error allocating arrays'

      inum = 0
      vol_conn(:,:) = -1
      do ispec = 1,NSPEC_CRUST_MANTLE
        do k = 1, NGLLZ-1
          do j = 1, NGLLY-1
            do i = 1, NGLLX-1
              ! indices of corner points
              id1 = vol_perm(ibool_crust_mantle(i,j,k,ispec))
              id2 = vol_perm(ibool_crust_mantle(i+1,j,k,ispec))
              id3 = vol_perm(ibool_crust_mantle(i+1,j+1,k,ispec))
              id4 = vol_perm(ibool_crust_mantle(i,j+1,k,ispec))

              id5 = vol_perm(ibool_crust_mantle(i,j,k+1,ispec))
              id6 = vol_perm(ibool_crust_mantle(i+1,j,k+1,ispec))
              id7 = vol_perm(ibool_crust_mantle(i+1,j+1,k+1,ispec))
              id8 = vol_perm(ibool_crust_mantle(i,j+1,k+1,ispec))

              ! note: indices for VTK start at 0
              inum = inum+1
              vol_conn(1,inum) = id1 - 1
              vol_conn(2,inum) = id2 - 1
              vol_conn(3,inum) = id3 - 1
              vol_conn(4,inum) = id4 - 1
              vol_conn(5,inum) = id5 - 1
              vol_conn(6,inum) = id6 - 1
              vol_conn(7,inum) = id7 - 1
              vol_conn(8,inum) = id8 - 1
            enddo
          enddo
        enddo
      enddo
    else
      ! point connectivity
      vol_nspec = NSPEC_CRUST_MANTLE

      allocate(vol_conn(8,vol_nspec),stat=ier)
      if (ier /= 0 ) stop 'Error allocating arrays'

      vol_conn(:,:) = -1
      do ispec = 1,NSPEC_CRUST_MANTLE
        ! indices of corner points
        id1 = vol_perm(ibool_crust_mantle(1,1,1,ispec))
        id2 = vol_perm(ibool_crust_mantle(NGLLX,1,1,ispec))
        id3 = vol_perm(ibool_crust_mantle(NGLLX,NGLLY,1,ispec))
        id4 = vol_perm(ibool_crust_mantle(1,NGLLY,1,ispec))

        id5 = vol_perm(ibool_crust_mantle(1,1,NGLLZ,ispec))
        id6 = vol_perm(ibool_crust_mantle(NGLLX,1,NGLLZ,ispec))
        id7 = vol_perm(ibool_crust_mantle(NGLLX,NGLLY,NGLLZ,ispec))
        id8 = vol_perm(ibool_crust_mantle(1,NGLLY,NGLLZ,ispec))

        ! note: indices for VTK start at 0
        vol_conn(1,ispec) = id1 - 1
        vol_conn(2,ispec) = id2 - 1
        vol_conn(3,ispec) = id3 - 1
        vol_conn(4,ispec) = id4 - 1
        vol_conn(5,ispec) = id5 - 1
        vol_conn(6,ispec) = id6 - 1
        vol_conn(7,ispec) = id7 - 1
        vol_conn(8,ispec) = id8 - 1
      enddo
    endif
    if (minval(vol_conn(:,:)) < 0) stop 'Error VTK volume point connectivity'

    ! allocates local data array
    allocate(vtkdata(vol_np),stat=ier)
    if (ier /= 0 ) stop 'Error allocating arrays'

    vtkdata(:) = 0.0

    ! gathers data from all MPI processes
    if (NPROC > 1) then
      ! multiple MPI processes

      ! user output
      !if (myrank == 0 ) print *,"    gathering all MPI info... "

      ! number of volume points for all partitions together
      call sum_all_i(vol_np,vtkdata_numpoints_all)
      if (myrank == 0 ) write(IMAIN,*) "    all volume points: ",vtkdata_numpoints_all

      ! gathers point info
      allocate(vtkdata_points_all(NPROC),stat=ier)
      if (ier /= 0 ) stop 'Error allocating arrays'

      vtkdata_points_all(:) = 0
      call gather_all_singlei(vol_np,vtkdata_points_all,NPROC)

      ! array offsets
      allocate(vtkdata_offset_all(NPROC),stat=ier)
      if (ier /= 0 ) stop 'Error allocating arrays'

      vtkdata_offset_all(1) = 0
      do i = 2, NPROC
        vtkdata_offset_all(i) = sum(vtkdata_points_all(1:i-1))
      enddo

      ! number of volume elements
      call sum_all_i(vol_nspec,vol_nspec_all)
      if (myrank == 0 ) write(IMAIN,*) "    all volume elements: ",vol_nspec_all

      ! volume elements
      allocate(vol_conn_nspec_all(NPROC),stat=ier)
      if (ier /= 0 ) stop 'Error allocating arrays'

      vol_conn_nspec_all(:) = 0
      call gather_all_singlei(8*vol_nspec,vol_conn_nspec_all,NPROC)

      ! array offsets
      allocate(vol_conn_offset_all(NPROC),stat=ier)
      if (ier /= 0 ) stop 'Error allocating arrays'

      vol_conn_offset_all(1) = 0
      do i = 2, NPROC
        vol_conn_offset_all(i) = sum(vol_conn_nspec_all(1:i-1))
      enddo

      ! global data arrays (only needed on main process)
      if (myrank == 0) then
        ! point data
        allocate(vtkdata_all(vtkdata_numpoints_all),stat=ier)
        if (ier /= 0 ) stop 'Error allocating vtkdata_all array'

        vtkdata_all(:) = 0.0

        ! gather locations
        allocate(vol_x_all(vtkdata_numpoints_all), &
                 vol_y_all(vtkdata_numpoints_all), &
                 vol_z_all(vtkdata_numpoints_all),stat=ier )
        if (ier /= 0 ) stop 'Error allocating vol_x_all,... arrays'

        vol_x_all(:) = 0.0
        vol_y_all(:) = 0.0
        vol_z_all(:) = 0.0

        ! connectivity
        allocate(vol_conn_all(8,vol_nspec_all),stat=ier)
        if (ier /= 0 ) stop 'Error allocating vol_conn_all array'

        vol_conn_all(:,:) = 0

      endif

      if (myrank == 0) then
        ! locations
        !if (myrank == 0 ) print *,"    locations..."
        call gatherv_all_r(vol_x,vol_np, &
                            vol_x_all,vtkdata_points_all,vtkdata_offset_all, &
                            vtkdata_numpoints_all,NPROC)
        call gatherv_all_r(vol_y,vol_np, &
                            vol_y_all,vtkdata_points_all,vtkdata_offset_all, &
                            vtkdata_numpoints_all,NPROC)
        call gatherv_all_r(vol_z,vol_np, &
                            vol_z_all,vtkdata_points_all,vtkdata_offset_all, &
                            vtkdata_numpoints_all,NPROC)

        ! connectivity
        !if (myrank == 0 ) print *,"    connectivity..."
        call gatherv_all_i(vol_conn,8*vol_nspec, &
                           vol_conn_all,vol_conn_nspec_all,vol_conn_offset_all, &
                           vol_nspec_all,NPROC)

        ! shifts connectivity ids for all additional slices
        do i = 2, NPROC
          ! divides by 8 to get nspec numbers
          ispec_start = vol_conn_offset_all(i)/8 + 1
          ispec_end = vol_conn_offset_all(i)/8 + vol_conn_nspec_all(i)/8
          do ispec = ispec_start,ispec_end
            vol_conn_all(:,ispec) = vol_conn_all(:,ispec) + vtkdata_offset_all(i)
          enddo
        enddo

        !if (myrank == 0 ) print *,"    preparing VTK field..."

        ! adds total volume wavefield to VTK window
        call prepare_vtkfield(vtkdata_numpoints_all,vol_x_all,vol_y_all,vol_z_all, &
                              vol_nspec_all,vol_conn_all)

      else
        ! all other process just send data
        ! locations
        call gatherv_all_r(vol_x,vol_np, &
                            dummy,vtkdata_points_all,vtkdata_offset_all, &
                            1,NPROC)
        call gatherv_all_r(vol_y,vol_np, &
                            dummy,vtkdata_points_all,vtkdata_offset_all, &
                            1,NPROC)
        call gatherv_all_r(vol_z,vol_np, &
                            dummy,vtkdata_points_all,vtkdata_offset_all, &
                            1,NPROC)
        ! connectivity
        call gatherv_all_i(vol_conn,8*vol_nspec, &
                            dummy_i,vol_conn_nspec_all,vol_conn_offset_all, &
                            1,NPROC)
      endif

    else
      ! serial run
      !if (myrank == 0 ) print *,"    preparing VTK field..."

      ! adds volume wavefield to VTK window
      call prepare_vtkfield(vol_np,vol_x,vol_y,vol_z,vol_nspec,vol_conn)
    endif

    ! frees memory
    deallocate(vol_x,vol_y,vol_z)
    deallocate(vol_conn,vol_perm)
    if (NPROC > 1) then
      deallocate(vol_conn_nspec_all,vol_conn_offset_all)
      if (myrank == 0 ) deallocate(vol_x_all,vol_y_all,vol_z_all,vol_conn_all)
    endif
  endif ! VTK_SHOW_VOLUME

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)"  VTK visualization preparation done"
    call flush_IMAIN()
  endif
  call synchronize_all()

  end subroutine prepare_vtk_window

