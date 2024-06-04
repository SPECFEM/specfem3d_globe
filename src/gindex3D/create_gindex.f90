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

  subroutine create_gindex()

  use gindex_par

  implicit none

  ! local parameters
  integer :: i_proc

  ! initialize global indices
  ignode_end = 0 ! global nodes for NGLLX = 5
  gnf_end = 0    ! global gdof for NGLLX = 5
  gnf_end1 = 0   ! global gdof for NGLLX_INF = 3

  !debug
  print *,'Number of solver processes: ',NPROCTOT

  ! folder for temporary files created by gindex3D
   if (myrank == 0) then
    write(IMAIN,*) 'creating temporary directory: tmp_gindex3D/'
    write(IMAIN,*)
    call flush_IMAIN()
  endif
  call execute_command_line('mkdir -p tmp_gindex3D/')

  ! loop through the processors
  do i_proc = 0,NPROCTOT-1
    ! creates global DOFs for this process
    call create_gindex_for_process(i_proc)
  enddo

  !debug
  print *,'all done'

  ! closes the main output file
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'End of gindex3D - all done'
    write(IMAIN,*)
    call flush_IMAIN()
    close(IMAIN)
  endif

  ! synchronizes all the processes to make sure everybody has finished
  call synchronize_all()

  end subroutine create_gindex

!
!------------------------------------------------------------------------
!

  subroutine create_gindex_for_process(i_proc)

  use gindex_par

  implicit none

  integer, intent(in) :: i_proc

  ! local parameters
  integer :: j_proc
  integer :: i,j,k,i_elmt,i_node,ier
  integer :: ispec_ic,ispec_oc,ispec_cm,ispec_trinf,ispec_inf
  integer :: ibool_ic,ibool_oc,ibool_cm,ibool_trinf,ibool_inf
  integer :: k_ic,k_oc,k_cm,k_trinf,k_inf
  integer :: ibool,inode,ignode,ispec,nnode_icb,nnode_cmb,nnode_trinfb,nnode_infb

  ! local
  integer,dimension(:),allocatable :: inode_ic,inode_oc,inode_cm,inode_trinf,inode_inf
  logical,dimension(:),allocatable :: isnode_ic,isnode_oc,isnode_cm,isnode_trinf,isnode_inf

  ! global
  integer,dimension(:),allocatable :: ignode_ic,ignode_oc,ignode_cm,ignode_trinf,ignode_inf
  logical,dimension(:),allocatable :: isgnode_ic,isgnode_oc,isgnode_cm,isgnode_trinf,isgnode_inf

  integer,allocatable :: inode_oc1(:),inode_ic1(:),inode_cm1(:),inode_trinf1(:),inode_inf1(:)
  integer,allocatable :: gdf_ic(:,:),gdf_oc(:,:),gdf_cm(:,:),gdf_trinf(:,:),gdf_inf(:,:),gnf(:,:)

  integer :: nnode_ic,nnode_oc,nnode_cm,nnode_trinf,nnode_inf
  integer :: inode1,inum,igll
  logical,allocatable :: isnode(:)

  logical,allocatable :: isibool_interface_ic(:,:),isibool_interface_oc(:,:), &
                         isibool_interface_cm(:,:),isibool_interface_trinf(:,:), &
                         isibool_interface_inf(:,:),isgnf(:,:)

  integer :: igdof,nibool
  integer,allocatable :: gghost(:,:),ighost(:)
  character(len = 20) :: fhead
  character(len = 12) :: spm,spn
  character(len = 128) :: fname

  integer,allocatable :: nmir(:)
  integer,allocatable :: gdf_ic1(:,:),gdf_oc1(:,:),gdf_cm1(:,:),gdf_trinf1(:,:), &
                         gdf_inf1(:,:),gnf1(:,:)

  logical,allocatable :: isgnf1(:,:)

  integer :: igdof1

  integer,allocatable :: tmpvec(:)
  integer,allocatable :: tmpmat(:,:)

  integer :: myrank_org

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'Process: ',i_proc
    write(IMAIN,*) '  reading databases...'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! to read databases for different processes, we will set myrank to the corresponding process id (i_proc)
  ! store current myrank
  myrank_org = myrank

  ! set new process id
  myrank = i_proc

  !debug
  print *,'Process: ',i_proc

  ! starts reading the databases
  call read_mesh_databases()

  ! restores original myrank
  myrank = myrank_org

  !deallocate unnecessary arrays
  ! inner core
  if (allocated(rmassz_inner_core)) then
    deallocate(rmassz_inner_core)
    if (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL) then
      deallocate(rmassx_inner_core,rmassy_inner_core)
    else
      nullify(rmassx_inner_core,rmassy_inner_core)
    endif
    deallocate(xstore_inner_core,ystore_inner_core,zstore_inner_core)
    deallocate(xix_inner_core,xiy_inner_core,xiz_inner_core, &
               etax_inner_core,etay_inner_core,etaz_inner_core, &
               gammax_inner_core,gammay_inner_core,gammaz_inner_core)
    deallocate(rhostore_inner_core,kappavstore_inner_core,muvstore_inner_core)
    deallocate(c11store_inner_core,c33store_inner_core,c12store_inner_core, &
               c13store_inner_core,c44store_inner_core)

    deallocate(phase_ispec_inner_inner_core)
    deallocate(num_elem_colors_inner_core)
    deallocate(buffer_send_vector_inner_core,buffer_recv_vector_inner_core, &
               request_send_vector_ic,request_recv_vector_ic)
  endif

  ! outer core
  if (allocated(rmass_outer_core)) then
    deallocate(rmass_outer_core)
    deallocate(xstore_outer_core,ystore_outer_core,zstore_outer_core)
    deallocate(xix_outer_core,xiy_outer_core,xiz_outer_core, &
               etax_outer_core,etay_outer_core,etaz_outer_core, &
               gammax_outer_core,gammay_outer_core,gammaz_outer_core)
    deallocate(rhostore_outer_core,kappavstore_outer_core)
    deallocate(vp_outer_core)

    deallocate(phase_ispec_inner_outer_core)
    deallocate(num_elem_colors_outer_core)
    deallocate(buffer_send_scalar_outer_core,buffer_recv_scalar_outer_core, &
               request_send_scalar_oc,request_recv_scalar_oc)
  endif

  ! crust/mantle
  if (allocated(rmassz_crust_mantle)) then
    deallocate(rmassz_crust_mantle)
    if ((NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS) .or. &
        (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL)) then
      deallocate(rmassx_crust_mantle,rmassy_crust_mantle)
    else
      nullify(rmassx_crust_mantle,rmassy_crust_mantle)
    endif
    deallocate(rmass_ocean_load)
    deallocate(xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle)
    deallocate(xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
               etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
               gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle)
    deallocate(rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle)
    deallocate(kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle)
    deallocate(c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle, &
               c14store_crust_mantle,c15store_crust_mantle,c16store_crust_mantle, &
               c22store_crust_mantle,c23store_crust_mantle,c24store_crust_mantle, &
               c25store_crust_mantle,c26store_crust_mantle,c33store_crust_mantle, &
               c34store_crust_mantle,c35store_crust_mantle,c36store_crust_mantle, &
               c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
               c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle)
    if (allocated(mu0store_crust_mantle)) deallocate(mu0store_crust_mantle)
    deallocate(ispec_is_tiso_crust_mantle)
    deallocate(rho_vp_crust_mantle,rho_vs_crust_mantle)

    deallocate(phase_ispec_inner_crust_mantle)
    deallocate(num_elem_colors_crust_mantle)
    deallocate(buffer_send_vector_crust_mantle,buffer_recv_vector_crust_mantle, &
               request_send_vector_cm,request_recv_vector_cm)
  endif

  ! trinfinite
  if (allocated(xstore_trinfinite)) then
    deallocate(xstore_trinfinite,ystore_trinfinite,zstore_trinfinite)
    deallocate(xix_trinfinite,xiy_trinfinite,xiz_trinfinite, &
               etax_trinfinite,etay_trinfinite,etaz_trinfinite, &
               gammax_trinfinite,gammay_trinfinite,gammaz_trinfinite)

    deallocate(phase_ispec_inner_trinfinite)
    deallocate(num_elem_colors_trinfinite)
    deallocate(buffer_send_scalar_trinfinite,buffer_recv_scalar_trinfinite, &
               request_send_scalar_trinfinite,request_recv_scalar_trinfinite)
  endif

  ! infinite
  if (allocated(xstore_infinite)) then
    deallocate(xstore_infinite,ystore_infinite,zstore_infinite)
    deallocate(xix_infinite,xiy_infinite,xiz_infinite, &
               etax_infinite,etay_infinite,etaz_infinite, &
               gammax_infinite,gammay_infinite,gammaz_infinite)

    deallocate(phase_ispec_inner_infinite)
    deallocate(num_elem_colors_infinite)
    deallocate(buffer_send_scalar_infinite,buffer_recv_scalar_infinite, &
               request_send_scalar_infinite,request_recv_scalar_infinite)
  endif

  ! coupling
  if (allocated(ibelm_moho_top)) then
    deallocate(ibelm_moho_top,ibelm_moho_bot, &
               ibelm_400_top,ibelm_400_bot, &
               ibelm_670_top,ibelm_670_bot, &
               normal_moho,normal_400,normal_670)
    deallocate(ibelm_xmin_crust_mantle,ibelm_xmax_crust_mantle, &
               ibelm_ymin_crust_mantle,ibelm_ymax_crust_mantle, &
               normal_xmin_crust_mantle,normal_xmax_crust_mantle, &
               normal_ymin_crust_mantle,normal_ymax_crust_mantle, &
               normal_bottom_crust_mantle,normal_top_crust_mantle, &
               jacobian2D_bottom_crust_mantle,jacobian2D_top_crust_mantle, &
               jacobian2D_xmin_crust_mantle,jacobian2D_xmax_crust_mantle, &
               jacobian2D_ymin_crust_mantle,jacobian2D_ymax_crust_mantle)
    deallocate(ibelm_xmin_outer_core,ibelm_xmax_outer_core, &
               ibelm_ymin_outer_core,ibelm_ymax_outer_core, &
               normal_xmin_outer_core,normal_xmax_outer_core, &
               normal_ymin_outer_core,normal_ymax_outer_core, &
               normal_bottom_outer_core,normal_top_outer_core, &
               jacobian2D_bottom_outer_core,jacobian2D_top_outer_core, &
               jacobian2D_xmin_outer_core,jacobian2D_xmax_outer_core, &
               jacobian2D_ymin_outer_core,jacobian2D_ymax_outer_core)
    deallocate(ibelm_xmin_inner_core,ibelm_xmax_inner_core, &
               ibelm_ymin_inner_core,ibelm_ymax_inner_core, &
               ibelm_bottom_inner_core)
  endif

  !debug
  !print *,'Process:',i_proc,' neighbors:',num_interfaces_inner_core

  nnode_ic = NGLOB_INNER_CORE
  nnode_oc = NGLOB_OUTER_CORE
  nnode_cm = NGLOB_CRUST_MANTLE
  nnode_trinf = NGLOB_TRINFINITE
  nnode_inf = NGLOB_INFINITE

  allocate(inode_ic(NGLOB_INNER_CORE), &
           inode_oc(NGLOB_OUTER_CORE), &
           inode_cm(NGLOB_CRUST_MANTLE), &
           inode_trinf(NGLOB_TRINFINITE), &
           inode_inf(NGLOB_INFINITE),stat=ier)
  if (ier /= 0) stop 'Error allocating inode_ic,.. arrays'
  inode_ic(:) = -1; inode_oc(:) = -1; inode_cm(:) = -1
  inode_trinf(:) = -1; inode_inf(:) = -1

  allocate(isnode_ic(NGLOB_INNER_CORE), &
           isnode_oc(NGLOB_OUTER_CORE), &
           isnode_cm(NGLOB_CRUST_MANTLE), &
           isnode_trinf(NGLOB_TRINFINITE), &
           isnode_inf(NGLOB_INFINITE),stat=ier)
  if (ier /= 0) stop 'Error allocating isnode_ic,.. arrays'
  isnode_ic(:) = .false.; isnode_oc(:) = .false.; isnode_cm(:) = .false.
  isnode_trinf(:) = .false.; isnode_inf(:) = .false.

  allocate(ignode_ic(NGLOB_INNER_CORE), &
           ignode_oc(NGLOB_OUTER_CORE), &
           ignode_cm(NGLOB_CRUST_MANTLE), &
           ignode_trinf(NGLOB_TRINFINITE), &
           ignode_inf(NGLOB_INFINITE),stat=ier)
  if (ier /= 0) stop 'Error allocating ignode_ic,.. arrays'
  ignode_ic(:) = -1; ignode_oc(:) = -1; ignode_cm(:) = -1
  ignode_trinf(:) = -1; ignode_inf(:) = -1

  allocate(isgnode_ic(NGLOB_INNER_CORE), &
           isgnode_oc(NGLOB_OUTER_CORE), &
           isgnode_cm(NGLOB_CRUST_MANTLE), &
           isgnode_trinf(NGLOB_TRINFINITE), &
           isgnode_inf(NGLOB_INFINITE),stat=ier)
  if (ier /= 0) stop 'Error allocating isgnode_ic,.. arrays'
  isgnode_ic(:) = .false.; isgnode_oc(:) = .false.; isgnode_cm(:) = .false.
  isgnode_trinf(:) = .false.; isgnode_inf(:) = .false.

  ! allocate necessary arrays
  if (.not. allocated(inode_elmt_cm)) then
    allocate(inode_elmt_cm(NGLLCUBE,NSPEC_CRUST_MANTLE))
    allocate(inode_elmt_cm1(NGLLCUBE_INF,NSPEC_CRUST_MANTLE))
  endif
  if (.not. allocated(inode_elmt_ic)) then
    allocate(inode_elmt_ic(NGLLCUBE,NSPEC_INNER_CORE))
    allocate(inode_elmt_ic1(NGLLCUBE_INF,NSPEC_INNER_CORE))
  endif
  if (.not. allocated(inode_elmt_oc)) then
    allocate(inode_elmt_oc(NGLLCUBE,NSPEC_OUTER_CORE))
    allocate(inode_elmt_oc1(NGLLCUBE_INF,NSPEC_OUTER_CORE))
  endif
  ! trinfinite arrays
  if (ADD_TRINF) then
    if (.not. allocated(inode_elmt_trinf)) then
      allocate(inode_elmt_trinf(NGLLCUBE,NSPEC_TRINFINITE))
      allocate(inode_elmt_trinf1(NGLLCUBE_INF,NSPEC_TRINFINITE))
    endif
  endif
  ! infinite arrays
  if (.not. allocated(inode_elmt_inf)) then
    allocate(inode_elmt_inf(NGLLCUBE,NSPEC_INFINITE))
    allocate(inode_elmt_inf1(NGLLCUBE_INF,NSPEC_INFINITE))
  endif

  ! count global node numbers
  nnode = NGLOB_INNER_CORE + NGLOB_OUTER_CORE + NGLOB_CRUST_MANTLE + NGLOB_TRINFINITE + NGLOB_INFINITE

  ! identify duplicate nodes on the boundary
  ! inner core - outer core boundary (ICB)
  !isnode_oc=.false.
  do i_elmt = 1,NSPEC2D_BOTTOM_OC
    ispec = ibelm_bottom_outer_core(i_elmt)
    k = 1 ! bottom face
    do j = 1,NGLLY
      do i = 1,NGLLX
        isgnode_oc(ibool_outer_core(i,j,k,ispec)) = .true.
      enddo
    enddo
  enddo
  nnode_icb = count(isgnode_oc)

  ! outer core - crust mantle boundary (CMB)
  !isnode_cm=.false.
  do i_elmt = 1,NSPEC2D_BOTTOM_CM
    ispec = ibelm_bottom_crust_mantle(i_elmt)
    k = 1 ! bottom face
    do j = 1,NGLLY
      do i = 1,NGLLX
        isgnode_cm(ibool_crust_mantle(i,j,k,ispec)) = .true.
      enddo
    enddo
  enddo
  nnode_cmb = count(isgnode_cm)

  ! crust mantle - transition infinite boundary (FS: free surface)
  nnode_trinfb = 0
  !isnode_trinf=.false.
  if (ADD_TRINF) then
    do i_elmt = 1,NSPEC2D_BOTTOM_TRINF
      ispec = ibelm_bottom_trinfinite(i_elmt)
      k = 1 ! bottom face
      do j = 1,NGLLY
        do i = 1,NGLLX
          isgnode_trinf(ibool_trinfinite(i,j,k,ispec)) = .true.
        enddo
      enddo
    enddo
    nnode_trinfb = count(isgnode_trinf)
  endif

  ! crust mantle - infinite boundary (FS: free surface)
  !isnode_inf=.false.
  do i_elmt = 1,NSPEC2D_BOTTOM_INF
    ispec = ibelm_bottom_infinite(i_elmt)
    k = 1 ! bottom face
    do j = 1,NGLLY
      do i = 1,NGLLX
        isgnode_inf(ibool_infinite(i,j,k,ispec)) = .true.
      enddo
    enddo
  enddo
  nnode_infb = count(isgnode_inf)

  ! number of unique global nodes in this processor
  nnode = nnode-nnode_icb-nnode_cmb-nnode_trinfb-nnode_infb

  print *,'Nodes in a process:',nnode

  ! until this stage both arrays are same
  !============================================================

  ! local indexing
  ! indexify regionally assembled local nodes and store in a region array
  ! inner core
  inode_ic = (/ (inode, inode = 1,NGLOB_INNER_CORE) /)

  inode = NGLOB_INNER_CORE
  ! outer core
  ! ICB
  ! copy common boundary nodes
  isnode_oc = .false.
  do i_elmt = 1,NSPEC2D_BOTTOM_OC
    ispec_oc = ibelm_bottom_outer_core(i_elmt)
    ispec_ic = ibelm_top_inner_core(i_elmt)
    k_oc = 1;    ! bottom face
    k_ic = NGLLZ ! top face
    do j = 1,NGLLY
      do i = 1,NGLLX
        ibool_oc = ibool_outer_core(i,j,k_oc,ispec_oc)
        ibool_ic = ibool_inner_core(i,j,k_ic,ispec_ic)
        inode_oc(ibool_oc)=inode_ic(ibool_ic)
        isnode_oc(ibool_oc) = .true.
      enddo
    enddo
  enddo
  ! now loop through all nodes
  do i_elmt = 1,NSPEC_OUTER_CORE
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          ibool_oc = ibool_outer_core(i,j,k,i_elmt)
           if (.not. isnode_oc(ibool_oc)) then
              isnode_oc(ibool_oc) = .true.
              inode = inode+1
              inode_oc(ibool_oc)=inode
           endif
        enddo
      enddo
    enddo
  enddo

  ! crust-mantle
  ! CMB
  ! copy common boundary nodes
  isnode_cm = .false.
  do i_elmt = 1,NSPEC2D_BOTTOM_CM
    ispec_cm = ibelm_bottom_crust_mantle(i_elmt)
    ispec_oc = ibelm_top_outer_core(i_elmt)
    k_cm = 1; k_oc = NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        ibool_cm = ibool_crust_mantle(i,j,k_cm,ispec_cm)
        ibool_oc = ibool_outer_core(i,j,k_oc,ispec_oc)
        inode_cm(ibool_cm)=inode_oc(ibool_oc)
        isnode_cm(ibool_cm) = .true.
      enddo
    enddo
  enddo
  ! now loop through all nodes
  do i_elmt = 1,NSPEC_CRUST_MANTLE
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          ibool_cm = ibool_crust_mantle(i,j,k,i_elmt)
          if (.not. isnode_cm(ibool_cm)) then
            isnode_cm(ibool_cm) = .true.
            inode = inode+1
            inode_cm(ibool_cm)=inode
          endif
        enddo
      enddo
    enddo
  enddo

  if (ADD_TRINF) then
    ! transition infinite & infinite
    ! FS
    ! copy common boundary nodes
    isnode_trinf = .false.
    do i_elmt = 1,NSPEC2D_BOTTOM_TRINF
      ispec_trinf = ibelm_bottom_trinfinite(i_elmt)
      ispec_cm = ibelm_top_crust_mantle(i_elmt)
      k_trinf = 1; k_cm = NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          ibool_trinf = ibool_trinfinite(i,j,k_trinf,ispec_trinf)
          ibool_cm = ibool_crust_mantle(i,j,k_cm,ispec_cm)
          inode_trinf(ibool_trinf)=inode_cm(ibool_cm)
          isnode_trinf(ibool_trinf) = .true.
        enddo
      enddo
    enddo
    ! now loop through all nodes
    do i_elmt = 1,NSPEC_TRINFINITE
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            ibool_trinf = ibool_trinfinite(i,j,k,i_elmt)
            if (.not. isnode_trinf(ibool_trinf)) then
              isnode_trinf(ibool_trinf) = .true.
              inode = inode+1
              inode_trinf(ibool_trinf)=inode
            endif
          enddo
        enddo
      enddo
    enddo

    ! infinite
    ! FS
    ! copy common boundary nodes
    isnode_inf = .false.
    do i_elmt = 1,NSPEC2D_BOTTOM_INF
      ispec_inf = ibelm_bottom_infinite(i_elmt)
      ispec_trinf = ibelm_top_trinfinite(i_elmt)
      k_inf = 1; k_trinf = NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          ibool_inf = ibool_infinite(i,j,k_inf,ispec_inf)
          ibool_trinf = ibool_trinfinite(i,j,k_trinf,ispec_trinf)
          inode_inf(ibool_inf)=inode_trinf(ibool_trinf)
          isnode_inf(ibool_inf) = .true.
        enddo
      enddo
    enddo
    ! now loop through all nodes
    do i_elmt = 1,NSPEC_INFINITE
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            ibool_inf = ibool_infinite(i,j,k,i_elmt)
            if (.not. isnode_inf(ibool_inf)) then
              isnode_inf(ibool_inf) = .true.
              inode = inode+1
              inode_inf(ibool_inf)=inode
            endif
          enddo
        enddo
      enddo
    enddo
  else
    ! infinite only
    ! FS
    ! copy common boundary nodes
    isnode_inf = .false.
    do i_elmt = 1,NSPEC2D_BOTTOM_INF
      ispec_inf = ibelm_bottom_infinite(i_elmt)
      ispec_cm = ibelm_top_crust_mantle(i_elmt)
      k_inf = 1; k_cm = NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          ibool_inf = ibool_infinite(i,j,k_inf,ispec_inf)
          ibool_cm = ibool_crust_mantle(i,j,k_cm,ispec_cm)
          inode_inf(ibool_inf)=inode_cm(ibool_cm)
          isnode_inf(ibool_inf) = .true.
        enddo
      enddo
    enddo
    ! now loop through all nodes
    do i_elmt = 1,NSPEC_INFINITE
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            ibool_inf = ibool_infinite(i,j,k,i_elmt)
            if (.not. isnode_inf(ibool_inf)) then
              isnode_inf(ibool_inf) = .true.
              inode = inode+1
              inode_inf(ibool_inf)=inode
            endif
          enddo
        enddo
      enddo
    enddo
  endif ! if (ADD_TRINF)

  ! safety check
  if (inode /= nnode) then
    write(*,*) 'ERROR: numbers of global nodes mismatch!',inode,nnode
    stop
  endif

  !============================================================

  ! global nodal indexing

  ! WARNING: is it correct to put these statements here?
  isgnode_ic(:) = .false.
  isgnode_oc(:) = .false.
  isgnode_cm(:) = .false.
  isgnode_trinf(:) = .false.
  isgnode_inf(:) = .false.

  ! copy global indices from preceding partitions

  ! inner core
  fhead = 'ic'
  write(spm,*) i_proc
  do i = 1,num_interfaces_inner_core
    j_proc = my_neighbors_inner_core(i)
    if (j_proc < i_proc) then
      write(spn,*) j_proc
      fname = 'tmp_gindex3D/'//trim(fhead)//trim(adjustl(spn))//'to'//trim(adjustl(spm))
      !print *,fname
      open(10,file=fname,action='read',status='old')
      read(10,*) nibool
      allocate(ighost(nibool))
      read(10,*) ighost(1:nibool)
      close(10,status='delete')

      isgnode_ic(ibool_interfaces_inner_core(1:nibool,i)) = .true.
      ignode_ic(ibool_interfaces_inner_core(1:nibool,i)) = ighost(:)
      deallocate(ighost)
    endif
  enddo

  ! outer core
  fhead = 'oc'
  write(spm,*) i_proc
  do i = 1,num_interfaces_outer_core
    j_proc = my_neighbors_outer_core(i)
    if (j_proc < i_proc) then
      write(spn,*) j_proc
      fname = 'tmp_gindex3D/'//trim(fhead)//trim(adjustl(spn))//'to'//trim(adjustl(spm))
      !print *,fname
      open(10,file=fname,action='read',status='old')
      read(10,*) nibool
      allocate(ighost(nibool))
      read(10,*) ighost(1:nibool)
      close(10,status='delete')

      isgnode_oc(ibool_interfaces_outer_core(1:nibool,i)) = .true.
      ignode_oc(ibool_interfaces_outer_core(1:nibool,i)) = ighost(:)
      deallocate(ighost)
    endif
  enddo

  ! crust mantle
  fhead = 'cm'
  write(spm,*) i_proc
  do i = 1,num_interfaces_crust_mantle
    j_proc = my_neighbors_crust_mantle(i)
    if (j_proc < i_proc) then
      write(spn,*) j_proc
      fname = 'tmp_gindex3D/'//trim(fhead)//trim(adjustl(spn))//'to'//trim(adjustl(spm))
      !print *,fname
      open(10,file=fname,action='read',status='old')
      read(10,*) nibool
      allocate(ighost(nibool))
      read(10,*) ighost(1:nibool)
      close(10,status='delete')

      isgnode_cm(ibool_interfaces_crust_mantle(1:nibool,i)) = .true.
      ignode_cm(ibool_interfaces_crust_mantle(1:nibool,i)) = ighost(:)
      deallocate(ighost)
    endif
  enddo

  if (ADD_TRINF) then
    ! transition infinite
    fhead = 'trinf'
    write(spm,*) i_proc
    do i = 1,num_interfaces_trinfinite
      j_proc = my_neighbors_trinfinite(i)
      if (j_proc < i_proc) then
        write(spn,*) j_proc
        fname = 'tmp_gindex3D/'//trim(fhead)//trim(adjustl(spn))//'to'//trim(adjustl(spm))
        !print *,fname
        open(10,file=fname,action='read',status='old')
        read(10,*) nibool
        allocate(ighost(nibool))
        read(10,*) ighost(1:nibool)
        close(10,status='delete')

        isgnode_trinf(ibool_interfaces_trinfinite(1:nibool,i)) = .true.
        ignode_trinf(ibool_interfaces_trinfinite(1:nibool,i)) = ighost(:)
        deallocate(ighost)
      endif
    enddo
  endif

  ! infinite
  fhead = 'inf'
  write(spm,*) i_proc
  do i = 1,num_interfaces_infinite
    j_proc = my_neighbors_infinite(i)
    if (j_proc < i_proc) then
      write(spn,*) j_proc
      fname = 'tmp_gindex3D/'//trim(fhead)//trim(adjustl(spn))//'to'//trim(adjustl(spm))
      !print *,fname
      open(10,file=fname,action='read',status='old')
      read(10,*) nibool
      allocate(ighost(nibool))
      read(10,*) ighost(1:nibool)
      close(10,status='delete')

      isgnode_inf(ibool_interfaces_infinite(1:nibool,i)) = .true.
      ignode_inf(ibool_interfaces_infinite(1:nibool,i)) = ighost(:)
      deallocate(ighost)
    endif
  enddo

  print *,'Previous largest node ID:',ignode_end

  ! indexify global nodes and store in a region array
  ! inner core
  ignode = ignode_end
  do i_node = 1,NGLOB_INNER_CORE
    if (.not. isgnode_ic(i_node)) then
      ignode = ignode+1
      isgnode_ic(i_node) = .true.
      ignode_ic(i_node) = ignode
    endif
  enddo
  if (i_proc == 1) print *,'IC:',ignode

  !ignode_ic=(/ (ignode, ignode=1,NGLOB_INNER_CORE) /)
  !ignode=NGLOB_INNER_CORE

  ! outer core
  ! ICB
  ! copy common boundary nodes
  !isnode_oc=.false.
  do i_elmt = 1,NSPEC2D_BOTTOM_OC
    ispec_oc = ibelm_bottom_outer_core(i_elmt)
    ispec_ic = ibelm_top_inner_core(i_elmt)
    k_oc = 1;    ! bottom face
    k_ic = NGLLZ ! top face
    do j = 1,NGLLY
      do i = 1,NGLLX
        ibool_oc = ibool_outer_core(i,j,k_oc,ispec_oc)
        ibool_ic = ibool_inner_core(i,j,k_ic,ispec_ic)
        ignode_oc(ibool_oc) = ignode_ic(ibool_ic)
        isgnode_oc(ibool_oc) = .true.
      enddo
    enddo
  enddo
  ! now loop through all nodes
  do i_elmt = 1,NSPEC_OUTER_CORE
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          ibool_oc = ibool_outer_core(i,j,k,i_elmt)
           if (.not. isgnode_oc(ibool_oc)) then
              isgnode_oc(ibool_oc) = .true.
              ignode = ignode+1
              ignode_oc(ibool_oc)=ignode
           endif
        enddo
      enddo
    enddo
  enddo
  if (i_proc == 1) print *,'OC:',ignode

  ! crust-mantle
  ! CMB
  ! copy common boundary nodes
  !isnode_cm=.false.
  do i_elmt = 1,NSPEC2D_BOTTOM_CM
    ispec_cm = ibelm_bottom_crust_mantle(i_elmt)
    ispec_oc = ibelm_top_outer_core(i_elmt)
    k_cm = 1; k_oc = NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        ibool_cm = ibool_crust_mantle(i,j,k_cm,ispec_cm)
        ibool_oc = ibool_outer_core(i,j,k_oc,ispec_oc)
        ignode_cm(ibool_cm) = ignode_oc(ibool_oc)
        isgnode_cm(ibool_cm) = .true.
      enddo
    enddo
  enddo
  ! now loop through all nodes
  do i_elmt = 1,NSPEC_CRUST_MANTLE
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          ibool_cm = ibool_crust_mantle(i,j,k,i_elmt)
          if (.not. isgnode_cm(ibool_cm)) then
            isgnode_cm(ibool_cm) = .true.
            ignode = ignode+1
            ignode_cm(ibool_cm)=ignode
          endif
        enddo
      enddo
    enddo
  enddo
  if (i_proc == 1) print *,'CM:',ignode

  if (ADD_TRINF) then
    ! transition infinite & infinite
    ! FS
    ! copy common boundary nodes
    !isnode_trinf=.false.
    do i_elmt = 1,NSPEC2D_BOTTOM_TRINF
      ispec_trinf = ibelm_bottom_trinfinite(i_elmt)
      ispec_cm = ibelm_top_crust_mantle(i_elmt)
      k_trinf = 1; k_cm = NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          ibool_trinf = ibool_trinfinite(i,j,k_trinf,ispec_trinf)
          ibool_cm = ibool_crust_mantle(i,j,k_cm,ispec_cm)
          ignode_trinf(ibool_trinf) = ignode_cm(ibool_cm)
          isgnode_trinf(ibool_trinf) = .true.
        enddo
      enddo
    enddo
    ! now loop through all nodes
    do i_elmt = 1,NSPEC_TRINFINITE
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            ibool_trinf = ibool_trinfinite(i,j,k,i_elmt)
            if (.not. isgnode_trinf(ibool_trinf)) then
              isgnode_trinf(ibool_trinf) = .true.
              ignode = ignode+1
              ignode_trinf(ibool_trinf) = ignode
            endif
          enddo
        enddo
      enddo
    enddo

    ! infinite
    ! FS
    ! copy common boundary nodes
    !isnode_inf=.false.
    do i_elmt = 1,NSPEC2D_BOTTOM_INF
      ispec_inf = ibelm_bottom_infinite(i_elmt)
      ispec_trinf = ibelm_top_trinfinite(i_elmt)
      k_inf = 1; k_trinf = NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          ibool_inf = ibool_infinite(i,j,k_inf,ispec_inf)
          ibool_trinf = ibool_trinfinite(i,j,k_trinf,ispec_trinf)
          ignode_inf(ibool_inf) = ignode_trinf(ibool_trinf)
          isgnode_inf(ibool_inf) = .true.
        enddo
      enddo
    enddo
    ! now loop through all nodes
    do i_elmt = 1,NSPEC_INFINITE
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            ibool_inf = ibool_infinite(i,j,k,i_elmt)
            if (.not. isgnode_inf(ibool_inf)) then
              isgnode_inf(ibool_inf) = .true.
              ignode = ignode+1
              ignode_inf(ibool_inf) = ignode
            endif
          enddo
        enddo
      enddo
    enddo
  else
    ! infinite only
    ! FS
    ! copy common boundary nodes
    !isnode_inf=.false.
    do i_elmt = 1,NSPEC2D_BOTTOM_INF
      ispec_inf = ibelm_bottom_infinite(i_elmt)
      ispec_cm = ibelm_top_crust_mantle(i_elmt)
      k_inf = 1; k_cm = NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          ibool_inf = ibool_infinite(i,j,k_inf,ispec_inf)
          ibool_cm = ibool_crust_mantle(i,j,k_cm,ispec_cm)
          ignode_inf(ibool_inf) = ignode_cm(ibool_cm)
          isgnode_inf(ibool_inf) = .true.
        enddo
      enddo
    enddo
    ! now loop through all nodes
    do i_elmt = 1,NSPEC_INFINITE
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            ibool_inf = ibool_infinite(i,j,k,i_elmt)
            if (.not. isgnode_inf(ibool_inf)) then
              isgnode_inf(ibool_inf) = .true.
              ignode = ignode+1
              ignode_inf(ibool_inf) = ignode
            endif
          enddo
        enddo
      enddo
    enddo
  endif
  if (i_proc == 1) print *,'INF:',ignode

  ! save global indices for neighboring partitions
  ! inner core
  fhead = 'ic'
  write(spm,*) i_proc
  do i = 1,num_interfaces_inner_core
    j_proc = my_neighbors_inner_core(i)
    if (j_proc > i_proc) then
      write(spn,*) j_proc
      fname = 'tmp_gindex3D/'//trim(fhead)//trim(adjustl(spm))//'to'//trim(adjustl(spn))
      !print *,fname !,i_proc,j_proc
      open(10,file=fname,action='write',status='replace')
      write(10,*) nibool_interfaces_inner_core(i)
      allocate(tmpvec(nibool_interfaces_inner_core(i)))
      tmpvec = ignode_ic(ibool_interfaces_inner_core(1:nibool_interfaces_inner_core(i),i))
      !write(10,*)ignode_ic(ibool_interfaces_inner_core(1:nibool_interfaces_inner_core(i),i))
      write(10,*) tmpvec
      deallocate(tmpvec)
      close(10)
    endif
  enddo

  ! outer core
  fhead = 'oc'
  write(spm,*) i_proc
  do i = 1,num_interfaces_outer_core
    j_proc = my_neighbors_outer_core(i)
    if (j_proc > i_proc) then
      write(spn,*) j_proc
      fname = 'tmp_gindex3D/'//trim(fhead)//trim(adjustl(spm))//'to'//trim(adjustl(spn))
      !print *,fname
      open(10,file=fname,action='write',status='replace')
      write(10,*) nibool_interfaces_outer_core(i)
      allocate(tmpvec(nibool_interfaces_outer_core(i)))
      tmpvec = ignode_oc(ibool_interfaces_outer_core(1:nibool_interfaces_outer_core(i),i))
      !write(10,*)ignode_oc(ibool_interfaces_outer_core(1:nibool_interfaces_outer_core(i),i))
      write(10,*) tmpvec
      deallocate(tmpvec)
      close(10)
    endif
  enddo

  ! crust mantle
  fhead = 'cm'
  write(spm,*) i_proc
  do i = 1,num_interfaces_crust_mantle
    j_proc = my_neighbors_crust_mantle(i)
    if (j_proc > i_proc) then
      write(spn,*) j_proc
      fname = 'tmp_gindex3D/'//trim(fhead)//trim(adjustl(spm))//'to'//trim(adjustl(spn))
      !print *,fname
      open(10,file=fname,action='write',status='replace')
      write(10,*) nibool_interfaces_crust_mantle(i)
      allocate(tmpvec(nibool_interfaces_crust_mantle(i)))
      tmpvec = ignode_cm(ibool_interfaces_crust_mantle(1:nibool_interfaces_crust_mantle(i),i))
      !write(10,*)ignode_cm(ibool_interfaces_crust_mantle(1:nibool_interfaces_crust_mantle(i),i))
      write(10,*) tmpvec
      deallocate(tmpvec)
      close(10)
    endif
  enddo

  if (ADD_TRINF) then
    ! transition infinite
    fhead = 'trinf'
    write(spm,*) i_proc
    do i = 1,num_interfaces_trinfinite
      j_proc = my_neighbors_trinfinite(i)
      if (j_proc > i_proc) then
        write(spn,*) j_proc
        fname = 'tmp_gindex3D/'//trim(fhead)//trim(adjustl(spm))//'to'//trim(adjustl(spn))
        !print *,fname
        open(10,file=fname,action='write',status='replace')
        write(10,*) nibool_interfaces_trinfinite(i)
        allocate(tmpvec(nibool_interfaces_trinfinite(i)))
        tmpvec = ignode_trinf(ibool_interfaces_trinfinite(1:nibool_interfaces_trinfinite(i),i))
        !write(10,*)ignode_trinf(ibool_interfaces_trinfinite(1:nibool_interfaces_trinfinite(i),i))
        write(10,*) tmpvec
        deallocate(tmpvec)
        close(10)
      endif
    enddo
  endif

  ! infinite
  fhead = 'inf'
  write(spm,*)i_proc
  do i = 1,num_interfaces_infinite
    j_proc = my_neighbors_infinite(i)
    if (j_proc > i_proc) then
      write(spn,*) j_proc
      fname = 'tmp_gindex3D/'//trim(fhead)//trim(adjustl(spm))//'to'//trim(adjustl(spn))
      !print *,fname
      open(10,file=fname,action='write',status='replace')
      write(10,*) nibool_interfaces_infinite(i)
      allocate(tmpvec(nibool_interfaces_infinite(i)))
      tmpvec = ignode_inf(ibool_interfaces_infinite(1:nibool_interfaces_infinite(i),i))
      !write(10,*)ignode_inf(ibool_interfaces_infinite(1:nibool_interfaces_infinite(i),i))
      write(10,*) tmpvec
      deallocate(tmpvec)
      close(10)
    endif
  enddo
  !  print *,'success!'

  ignode_end = maxval( (/ ignode_ic,ignode_oc,ignode_cm,ignode_trinf,ignode_inf /) )
  print *,'Largest node ID:',ignode_end

  write(spm,*) i_proc

  ! file output
  fname = 'DATABASES_MPI/gibool_proc'//trim(adjustl(spm))
  open(10,file=fname,action='write',status='replace')

  write(10,*) NGLOB_INNER_CORE
  write(10,*) ignode_ic
  write(10,*) NGLOB_OUTER_CORE
  write(10,*) ignode_oc
  write(10,*) NGLOB_CRUST_MANTLE
  write(10,*) ignode_cm
  write(10,*) NGLOB_TRINFINITE
  write(10,*) ignode_trinf
  write(10,*) NGLOB_INFINITE
  write(10,*) ignode_inf
  close(10)

  ! global indexing of degrees of freedoms
  allocate(gnf(NNDOF,nnode),isgnf(NNDOF,nnode))
  gnf = 0
  isgnf = .false.

  ! activate freedoms

  ! freedoms of fictitious cube in inner core are deactivated
  do i_elmt = 1,NSPEC_INNER_CORE
    if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE) cycle
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          ibool = ibool_inner_core(i,j,k,i_elmt)
          gnf(1,inode_ic(ibool)) = 1
        enddo
      enddo
    enddo
  enddo

  ! outer core
  ! all freedoms are active
  gnf(1,inode_oc) = 1

  ! crust-mantle
  ! all freedoms are active
  gnf(1,inode_cm) = 1

  ! transition infinite
  if (ADD_TRINF) then
    ! all freedoms are active
    gnf(1,inode_trinf) = 1
  endif

  ! infinite element
  ! all but surface nodes are activated
  do i_elmt = 1,NSPEC_INFINITE
    do k = 1,NGLLZ-1
      do j = 1,NGLLY
        do i = 1,NGLLX
          ibool = ibool_infinite(i,j,k,i_elmt)
          gnf(1,inode_inf(ibool)) = 1
        enddo
      enddo
    enddo
  enddo

  ! copy global indices from preceding partitions
  ! inner core
  fhead = 'gdof_ic'
  write(spm,*) i_proc
  do i = 1,num_interfaces_inner_core
    j_proc = my_neighbors_inner_core(i)
    if (j_proc < i_proc) then
      write(spn,*) j_proc
      fname = 'tmp_gindex3D/'//trim(fhead)//trim(adjustl(spn))//'to'//trim(adjustl(spm))
      !print *,fname
      open(10,file=fname,action='read',status='old')
      read(10,*) nibool
      allocate(gghost(NNDOF,nibool))
      read(10,*) gghost
      close(10,status='delete')

      isgnf(:,inode_ic(ibool_interfaces_inner_core(1:nibool,i))) = .true.
      gnf(:,inode_ic(ibool_interfaces_inner_core(1:nibool,i))) = gghost(:,:)
      deallocate(gghost)
    endif
  enddo

  ! outer core
  fhead = 'gdof_oc'
  write(spm,*) i_proc
  do i = 1,num_interfaces_outer_core
    j_proc = my_neighbors_outer_core(i)
    if (j_proc < i_proc) then
      write(spn,*) j_proc
      fname = 'tmp_gindex3D/'//trim(fhead)//trim(adjustl(spn))//'to'//trim(adjustl(spm))
      !print *,fname
      open(10,file=fname,action='read',status='old')
      read(10,*) nibool
      allocate(gghost(NNDOF,nibool))
      read(10,*) gghost
      close(10,status='delete')

      isgnf(:,inode_oc(ibool_interfaces_outer_core(1:nibool,i))) = .true.
      gnf(:,inode_oc(ibool_interfaces_outer_core(1:nibool,i))) = gghost(:,:)
      deallocate(gghost)
    endif
  enddo

  ! crust mantle
  fhead = 'gdof_cm'
  write(spm,*) i_proc
  do i = 1,num_interfaces_crust_mantle
    j_proc = my_neighbors_crust_mantle(i)
    if (j_proc < i_proc) then
      write(spn,*) j_proc
      fname = 'tmp_gindex3D/'//trim(fhead)//trim(adjustl(spn))//'to'//trim(adjustl(spm))
      !print *,fname
      open(10,file=fname,action='read',status='old')
      read(10,*) nibool
      allocate(gghost(NNDOF,nibool))
      read(10,*) gghost
      close(10,status='delete')

      isgnf(:,inode_cm(ibool_interfaces_crust_mantle(1:nibool,i))) = .true.
      gnf(:,inode_cm(ibool_interfaces_crust_mantle(1:nibool,i))) = gghost(:,:)
      deallocate(gghost)
    endif
  enddo

  if (ADD_TRINF) then
    ! transition infinite
    fhead = 'gdof_trinf'
    write(spm,*) i_proc
    do i = 1,num_interfaces_trinfinite
      j_proc = my_neighbors_trinfinite(i)
      if (j_proc < i_proc) then
        write(spn,*) j_proc
        fname = 'tmp_gindex3D/'//trim(fhead)//trim(adjustl(spn))//'to'//trim(adjustl(spm))
        !print *,fname
        open(10,file=fname,action='read',status='old')
        read(10,*) nibool
        allocate(gghost(NNDOF,nibool))
        read(10,*) gghost
        close(10,status='delete')

        isgnf(:,inode_trinf(ibool_interfaces_trinfinite(1:nibool,i))) = .true.
        gnf(:,inode_trinf(ibool_interfaces_trinfinite(1:nibool,i))) = gghost(:,:)
        deallocate(gghost)
      endif
    enddo
  endif

  ! infinite
  fhead = 'gdof_inf'
  write(spm,*) i_proc
  do i = 1,num_interfaces_infinite
    j_proc = my_neighbors_infinite(i)
    if (j_proc < i_proc) then
      write(spn,*) j_proc
      fname = 'tmp_gindex3D/'//trim(fhead)//trim(adjustl(spn))//'to'//trim(adjustl(spm))
      !print *,fname
      open(10,file=fname,action='read',status='old')
      read(10,*) nibool
      allocate(gghost(NNDOF,nibool))
      read(10,*) gghost
      close(10,status='delete')

      isgnf(:,inode_inf(ibool_interfaces_infinite(1:nibool,i))) = .true.
      gnf(:,inode_inf(ibool_interfaces_infinite(1:nibool,i))) = gghost(:,:)
      deallocate(gghost)
    endif
  enddo

  print *,'Previous largest gnf ID:',gnf_end

  igdof = gnf_end ! gdof
  do i_node = 1,nnode
    if (gnf(1,i_node) > 0 .and. .not. isgnf(1,i_node)) then
      isgnf(1,i_node) = .true.
      igdof = igdof+1
      gnf(1,i_node) = igdof
    endif
  enddo
  !neq=inode

  allocate(gdf_ic(NNDOF,NGLOB_INNER_CORE), &
           gdf_oc(NNDOF,NGLOB_OUTER_CORE), &
           gdf_cm(NNDOF,NGLOB_CRUST_MANTLE), &
           gdf_trinf(NNDOF,NGLOB_TRINFINITE), &
           gdf_inf(NNDOF,NGLOB_INFINITE))
  ! store gdf in a region array
  gdf_ic(:,:) = gnf(:,inode_ic)
  gdf_oc(:,:) = gnf(:,inode_oc)
  gdf_cm(:,:) = gnf(:,inode_cm)
  gdf_trinf(:,:) = gnf(:,inode_trinf)
  gdf_inf(:,:) = gnf(:,inode_inf)

  ! save global degrees of freedom for neighboring partitions
  ! inner core
  fhead = 'gdof_ic'
  write(spm,*) i_proc
  do i = 1,num_interfaces_inner_core
    j_proc = my_neighbors_inner_core(i)
    if (j_proc > i_proc) then
      write(spn,*) j_proc
      fname = 'tmp_gindex3D/'//trim(fhead)//trim(adjustl(spm))//'to'//trim(adjustl(spn))
      !print *,fname !,i_proc,j_proc
      open(10,file=fname,action='write',status='replace')
      write(10,*) nibool_interfaces_inner_core(i)
      allocate(tmpmat(NNDOF,nibool_interfaces_inner_core(i)))
      tmpmat = gdf_ic(:,ibool_interfaces_inner_core(1:nibool_interfaces_inner_core(i),i))
      !write(10,*)gdf_ic(:,ibool_interfaces_inner_core(1:nibool_interfaces_inner_core(i),i))
      write(10,*) tmpmat
      deallocate(tmpmat)
      close(10)
    endif
  enddo

  ! outer core
  fhead = 'gdof_oc'
  write(spm,*) i_proc
  do i = 1,num_interfaces_outer_core
    j_proc = my_neighbors_outer_core(i)
    if (j_proc > i_proc) then
      write(spn,*) j_proc
      fname = 'tmp_gindex3D/'//trim(fhead)//trim(adjustl(spm))//'to'//trim(adjustl(spn))
      !print *,fname
      open(10,file=fname,action='write',status='replace')
      write(10,*)nibool_interfaces_outer_core(i)
      allocate(tmpmat(NNDOF,nibool_interfaces_outer_core(i)))
      tmpmat = gdf_oc(:,ibool_interfaces_outer_core(1:nibool_interfaces_outer_core(i),i))
      !write(10,*)gdf_oc(:,ibool_interfaces_outer_core(1:nibool_interfaces_outer_core(i),i))
      write(10,*) tmpmat
      deallocate(tmpmat)
      close(10)
    endif
  enddo

  ! crust mantle
  fhead = 'gdof_cm'
  write(spm,*) i_proc
  do i = 1,num_interfaces_crust_mantle
    j_proc = my_neighbors_crust_mantle(i)
    if (j_proc > i_proc) then
      write(spn,*) j_proc
      fname = 'tmp_gindex3D/'//trim(fhead)//trim(adjustl(spm))//'to'//trim(adjustl(spn))
      !print *,fname
      open(10,file=fname,action='write',status='replace')
      write(10,*) nibool_interfaces_crust_mantle(i)
      allocate(tmpmat(NNDOF,nibool_interfaces_crust_mantle(i)))
      tmpmat = gdf_cm(:,ibool_interfaces_crust_mantle(1:nibool_interfaces_crust_mantle(i),i))
      !write(10,*)gdf_cm(:,ibool_interfaces_crust_mantle(1:nibool_interfaces_crust_mantle(i),i))
      write(10,*) tmpmat
      deallocate(tmpmat)
      close(10)
    endif
  enddo

  if (ADD_TRINF) then
    ! transition infinite
    fhead = 'gdof_trinf'
    write(spm,*) i_proc
    do i = 1,num_interfaces_trinfinite
      j_proc = my_neighbors_trinfinite(i)
      if (j_proc > i_proc) then
        write(spn,*) j_proc
        fname = 'tmp_gindex3D/'//trim(fhead)//trim(adjustl(spm))//'to'//trim(adjustl(spn))
        !print *,fname
        open(10,file=fname,action='write',status='replace')
        write(10,*) nibool_interfaces_trinfinite(i)
        allocate(tmpmat(NNDOF,nibool_interfaces_trinfinite(i)))
        tmpmat = gdf_trinf(:,ibool_interfaces_trinfinite(1:nibool_interfaces_trinfinite(i),i))
        !write(10,*)gdf_trinf(:,ibool_interfaces_trinfinite(1:nibool_interfaces_trinfinite(i),i))
        write(10,*) tmpmat
        deallocate(tmpmat)
        close(10)
      endif
    enddo
  endif

  ! infinite
  fhead = 'gdof_inf'
  write(spm,*)i_proc
  do i = 1,num_interfaces_infinite
    j_proc = my_neighbors_infinite(i)
    if (j_proc > i_proc) then
      write(spn,*) j_proc
      fname = 'tmp_gindex3D/'//trim(fhead)//trim(adjustl(spm))//'to'//trim(adjustl(spn))
      !print *,fname
      open(10,file=fname,action='write',status='replace')
      write(10,*) nibool_interfaces_infinite(i)
      allocate(tmpmat(NNDOF,nibool_interfaces_infinite(i)))
      tmpmat = gdf_inf(:,ibool_interfaces_infinite(1:nibool_interfaces_infinite(i),i))
      !write(10,*)gdf_inf(:,ibool_interfaces_infinite(1:nibool_interfaces_infinite(i),i))
      write(10,*) tmpmat
      deallocate(tmpmat)
      close(10)
    endif
  enddo

  gnf_end = maxval(gnf)
  print *,'Largest gnf ID:',gnf_end

  write(spm,*) i_proc

  ! file output
  ! (needed for PETSc Level-2 solver setup)
  fname='DATABASES_MPI/gdof_proc'//trim(adjustl(spm))
  open(10,file=fname,action='write',status='replace')

  write(10,*) NGLOB_INNER_CORE
  write(10,*) gdf_ic
  write(10,*) NGLOB_OUTER_CORE
  write(10,*) gdf_oc
  write(10,*) NGLOB_CRUST_MANTLE
  write(10,*) gdf_cm
  write(10,*) NGLOB_TRINFINITE
  write(10,*) gdf_trinf
  write(10,*) NGLOB_INFINITE
  write(10,*) gdf_inf
  close(10)

  !=============================================================================
  ! global degrees of freedoms for NGLLX_INF=3
  !=============================================================================
  !
  ! activate GLL points for NGLLX_INF=3 from NGLLX=5
  is_active_gll(:) = .false.
  inum = 0
  do k = 1,NGLLZ,2
    do j = 1,NGLLY,2
      do i = 1,NGLLX,2
        inum = inum+1
        igll = NGLLY * NGLLX * (k-1) + NGLLX * (j-1) + i
        is_active_gll(igll) = .true.
        igll_active_on(inum) = igll
      enddo
    enddo
  enddo

  ! prepare mesh information
  ! store ibool in 1D linear mapping
  ! inner core
  do i_elmt = 1,NSPEC_INNER_CORE
    inode_elmt_ic(:,i_elmt) = reshape(ibool_inner_core(:,:,:,i_elmt),(/NGLLCUBE/))
  enddo
  ! outer core
  do i_elmt = 1,NSPEC_OUTER_CORE
    inode_elmt_oc(:,i_elmt) = reshape(ibool_outer_core(:,:,:,i_elmt),(/NGLLCUBE/))
  enddo
  ! inner core
  do i_elmt = 1,NSPEC_CRUST_MANTLE
    inode_elmt_cm(:,i_elmt) = reshape(ibool_crust_mantle(:,:,:,i_elmt),(/NGLLCUBE/))
  enddo
  ! transition infinite
  if (ADD_TRINF) then
  do i_elmt = 1,NSPEC_TRINFINITE
    inode_elmt_trinf(:,i_elmt) = reshape(ibool_trinfinite(:,:,:,i_elmt),(/NGLLCUBE/))
  enddo
  endif
  ! infinite
  do i_elmt = 1,NSPEC_INFINITE
    inode_elmt_inf(:,i_elmt) = reshape(ibool_infinite(:,:,:,i_elmt),(/NGLLCUBE/))
  enddo

  isnode_ic(:) = .false.
  isnode_oc(:) = .false.
  isnode_cm(:) = .false.
  isnode_trinf(:) = .false.
  isnode_inf(:) = .false.

  allocate(isnode(nnode))
  isnode(:) = .false.

  ! inner core
  do i_elmt = 1,NSPEC_INNER_CORE
    do k = 1,NGLLZ,2
      do j = 1,NGLLY,2
        do i = 1,NGLLX,2
          ibool = ibool_inner_core(i,j,k,i_elmt)
          isnode_ic(ibool) = .true.
          isnode(inode_ic(ibool)) = .true.
        enddo
      enddo
    enddo
  enddo

  ! outer core
  do i_elmt = 1,NSPEC_OUTER_CORE
    do k = 1,NGLLZ,2
      do j = 1,NGLLY,2
        do i = 1,NGLLX,2
          ibool = ibool_outer_core(i,j,k,i_elmt)
          isnode_oc(ibool) = .true.
          isnode(inode_oc(ibool)) = .true.
        enddo
      enddo
    enddo
  enddo

  ! crust mantle
  do i_elmt = 1,NSPEC_CRUST_MANTLE
    do k = 1,NGLLZ,2
      do j = 1,NGLLY,2
        do i = 1,NGLLX,2
          ibool = ibool_crust_mantle(i,j,k,i_elmt)
          isnode_cm(ibool) = .true.
          isnode(inode_cm(ibool)) = .true.
        enddo
      enddo
    enddo
  enddo

  if (ADD_TRINF) then
    ! trinfinite
    do i_elmt = 1,NSPEC_TRINFINITE
      do k = 1,NGLLZ,2
        do j = 1,NGLLY,2
          do i = 1,NGLLX,2
            ibool = ibool_trinfinite(i,j,k,i_elmt)
            isnode_trinf(ibool) = .true.
            isnode(inode_trinf(ibool)) = .true.
          enddo
        enddo
      enddo
    enddo
  endif

  ! infinite
  do i_elmt = 1,NSPEC_INFINITE
    do k = 1,NGLLZ,2
      do j = 1,NGLLY,2
        do i = 1,NGLLX,2
          ibool = ibool_infinite(i,j,k,i_elmt)
          isnode_inf(ibool) = .true.
          isnode(inode_inf(ibool)) = .true.
        enddo
      enddo
    enddo
  enddo

  nnode_ic1 = count(isnode_ic)
  nnode_oc1 = count(isnode_oc)
  nnode_cm1 = count(isnode_cm)
  nnode_trinf1 = count(isnode_trinf)
  nnode_inf1 = count(isnode_inf)

  nnode1 = count(isnode)

  ! node mirror
  allocate(nmir(nnode))
  inode1 = 0
  do i_node = 1,nnode
    if (isnode(i_node)) then
      inode1 = inode1+1
      nmir(i_node) = inode1
    endif
  enddo
  deallocate(isnode)

  if (inode1 /= nnode1) then
    print *,inode1,nnode1
    write(*,'(/,a) ') 'ERROR: counted level-1 active nodes mismatch!'
    stop
  endif

  ! find active nodes and mirror to orginal nodes
  ! inner core
  allocate(nmir_ic(nnode_ic))
  nmir_ic(:) = 0
  inode = 0
  do i_node = 1,nnode_ic
    if (isnode_ic(i_node)) then
      inode = inode+1
      nmir_ic(i_node) = inode
    endif
  enddo
  if (inode /= nnode_ic1) then
    print *,inode,nnode_ic1,nnode_ic,size(isnode_ic),NGLOB_INNER_CORE
    write(*,'(/,a) ') 'ERROR: counted level-1 active nodes mismatch in inner core!'
    stop
  endif

  ! outer core
  allocate(nmir_oc(nnode_oc))
  nmir_oc(:) = 0
  inode = 0
  do i_node = 1,nnode_oc
    if (isnode_oc(i_node)) then
      inode = inode+1
      nmir_oc(i_node) = inode
    endif
  enddo
  if (inode /= nnode_oc1) then
    write(*,'(/,a) ') 'ERROR: counted level-1 active nodes mismatch in outer core!'
    stop
  endif

  ! crust mantle
  allocate(nmir_cm(nnode_cm))
  nmir_cm(:) = 0
  inode = 0
  do i_node = 1,nnode_cm
    if (isnode_cm(i_node)) then
      inode = inode+1
      nmir_cm(i_node) = inode
    endif
  enddo
  if (inode /= nnode_cm1) then
    write(*,'(/,a) ') 'ERROR: counted level-1 active nodes mismatch in crust mantle!'
    stop
  endif

  ! transition infinite
  if (ADD_TRINF) then
    allocate(nmir_trinf(nnode_trinf))
    nmir_trinf(:) = 0
    inode = 0
    do i_node = 1,nnode_trinf
      if (isnode_trinf(i_node)) then
        inode = inode+1
        nmir_trinf(i_node) = inode
      endif
    enddo
    if (inode /= nnode_trinf1) then
      write(*,'(/,a) ') 'ERROR: counted level-1 active nodes mismatch in transition infinite!'
      stop
    endif
  endif

  ! infinite
  allocate(nmir_inf(nnode_inf))
  nmir_inf(:) = 0
  inode = 0
  do i_node = 1,nnode_inf
    if (isnode_inf(i_node)) then
      inode = inode+1
      nmir_inf(i_node) = inode
    endif
  enddo
  if (inode /= nnode_inf1) then
    write(*,'(/,a) ') 'ERROR: counted level-1 active nodes mismatch in infinite!'
    stop
  endif

  ! store ibool1 in elemental array
  ! inner core
  do i_elmt = 1,NSPEC_INNER_CORE
    inode_elmt_ic1(:,i_elmt) = nmir_ic(inode_elmt_ic(igll_active_on,i_elmt))
  enddo

  ! outer core
  do i_elmt = 1,NSPEC_OUTER_CORE
    inode_elmt_oc1(:,i_elmt) = nmir_oc(inode_elmt_oc(igll_active_on,i_elmt))
  enddo

  ! crust mantle
  do i_elmt = 1,NSPEC_CRUST_MANTLE
    inode_elmt_cm1(:,i_elmt) = nmir_cm(inode_elmt_cm(igll_active_on,i_elmt))
  enddo

  ! transtion infinite
  if (ADD_TRINF) then
    do i_elmt = 1,NSPEC_TRINFINITE
      inode_elmt_trinf1(:,i_elmt) = nmir_trinf(inode_elmt_trinf(igll_active_on,i_elmt))
    enddo
  endif

  ! infinite
  do i_elmt = 1,NSPEC_INFINITE
    inode_elmt_inf1(:,i_elmt) = nmir_inf(inode_elmt_inf(igll_active_on,i_elmt))
  enddo

  ! find inode_ic1,inode_oc1,inode_cm1,inode_inf1
  allocate(inode_ic1(nnode_ic1),inode_oc1(nnode_oc1),inode_cm1(nnode_cm1), &
  inode_trinf1(nnode_trinf1),inode_inf1(nnode_inf1))

  ! inner core
  do i_node = 1,nnode_ic
    if (isnode_ic(i_node)) then
      inode_ic1(nmir_ic(i_node)) = nmir(inode_ic(i_node))
    endif
  enddo

  ! outer core
  do i_node = 1,nnode_oc
    if (isnode_oc(i_node)) then
      inode_oc1(nmir_oc(i_node)) = nmir(inode_oc(i_node))
    endif
  enddo

  ! crust mantle
  do i_node = 1,nnode_cm
    if (isnode_cm(i_node)) then
      inode_cm1(nmir_cm(i_node)) = nmir(inode_cm(i_node))
    endif
  enddo

  ! transition infinite
  if (ADD_TRINF) then
    do i_node = 1,nnode_trinf
      if (isnode_trinf(i_node)) then
        inode_trinf1(nmir_trinf(i_node)) = nmir(inode_trinf(i_node))
      endif
    enddo
  endif

  ! infinite
  do i_node = 1,nnode_inf
    if (isnode_inf(i_node)) then
      inode_inf1(nmir_inf(i_node))=nmir(inode_inf(i_node))
    endif
  enddo
  deallocate(nmir)

  ! modify MPI interfaces removing deactivated nodes
  ! NOTE: this has to be done before isnode_ic, isnode_oc, etc. change

  ! inner core
  num_interfaces_inner_core1 = num_interfaces_inner_core
  allocate(my_neighbors_inner_core1(num_interfaces_inner_core1))
  if (num_interfaces_inner_core1 > 0) my_neighbors_inner_core1(:) = my_neighbors_inner_core(:)

  allocate(isibool_interface_ic(max_nibool_interfaces_ic,num_interfaces_inner_core1), &
           nibool_interfaces_inner_core1(num_interfaces_inner_core1))
  isibool_interface_ic(:,:) = .false.
  nibool_interfaces_inner_core1(:) = 0
  do i = 1,num_interfaces_inner_core1
    do j = 1,nibool_interfaces_inner_core(i)
      isibool_interface_ic(j,i) = isnode_ic(ibool_interfaces_inner_core(j,i))
    enddo
    nibool_interfaces_inner_core1(i) = count(isibool_interface_ic(:,i))
  enddo
  if (num_interfaces_inner_core1 > 0) then
    max_nibool_interfaces_inner_core1 = maxval(nibool_interfaces_inner_core1)
  else
    max_nibool_interfaces_inner_core1 = 0
  endif

  allocate(ibool_interfaces_inner_core1(max_nibool_interfaces_inner_core1,num_interfaces_inner_core1))
  ibool_interfaces_inner_core1(:,:) = 0
  do i = 1,num_interfaces_inner_core1
    inum = 0
    do j = 1,nibool_interfaces_inner_core(i)
      if (isibool_interface_ic(j,i)) then
        inum = inum+1
        ibool_interfaces_inner_core1(inum,i) = nmir_ic(ibool_interfaces_inner_core(j,i))
      endif
    enddo
  enddo

  deallocate(nmir_ic)

  ! outer core
  num_interfaces_outer_core1 = num_interfaces_outer_core
  allocate(my_neighbors_outer_core1(num_interfaces_outer_core1))
  if (num_interfaces_outer_core1 > 0) my_neighbors_outer_core1(:) = my_neighbors_outer_core(:)

  allocate(isibool_interface_oc(max_nibool_interfaces_oc,num_interfaces_outer_core1), &
           nibool_interfaces_outer_core1(num_interfaces_outer_core1))
  isibool_interface_oc(:,:) = .false.
  nibool_interfaces_outer_core1(:) = 0
  do i = 1,num_interfaces_outer_core1
    do j = 1,nibool_interfaces_outer_core(i)
      isibool_interface_oc(j,i)=isnode_oc(ibool_interfaces_outer_core(j,i))
    enddo
    nibool_interfaces_outer_core1(i)=count(isibool_interface_oc(:,i))
  enddo
  if (num_interfaces_outer_core1 > 0) then
    max_nibool_interfaces_outer_core1 = maxval(nibool_interfaces_outer_core1)
  else
    max_nibool_interfaces_outer_core1 = 0
  endif

  allocate(ibool_interfaces_outer_core1(max_nibool_interfaces_outer_core1,num_interfaces_outer_core1))
  ibool_interfaces_outer_core1(:,:) = 0
  do i = 1,num_interfaces_outer_core1
    inum = 0
    do j = 1,nibool_interfaces_outer_core(i)
      if (isibool_interface_oc(j,i)) then
        inum = inum+1
        ibool_interfaces_outer_core1(inum,i) = nmir_oc(ibool_interfaces_outer_core(j,i))
      endif
    enddo
  enddo
  deallocate(nmir_oc)

  ! crust mantle
  num_interfaces_crust_mantle1 = num_interfaces_crust_mantle
  allocate(my_neighbors_crust_mantle1(num_interfaces_crust_mantle1))
  if (num_interfaces_crust_mantle1 > 0) my_neighbors_crust_mantle1(:) = my_neighbors_crust_mantle(:)

  allocate(isibool_interface_cm(max_nibool_interfaces_cm,num_interfaces_crust_mantle1), &
           nibool_interfaces_crust_mantle1(num_interfaces_crust_mantle1))
  isibool_interface_cm(:,:) = .false.
  nibool_interfaces_crust_mantle1(:) = 0
  do i = 1,num_interfaces_crust_mantle1
    do j = 1,nibool_interfaces_crust_mantle(i)
      isibool_interface_cm(j,i)=isnode_cm(ibool_interfaces_crust_mantle(j,i))
    enddo
    nibool_interfaces_crust_mantle1(i)=count(isibool_interface_cm(:,i))
  enddo
  if (num_interfaces_crust_mantle1 > 0) then
    max_nibool_interfaces_crust_mantle1 = maxval(nibool_interfaces_crust_mantle1)
  else
    max_nibool_interfaces_crust_mantle1 = 0
  endif

  allocate(ibool_interfaces_crust_mantle1(max_nibool_interfaces_crust_mantle1,num_interfaces_crust_mantle1))
  ibool_interfaces_crust_mantle1(:,:) = 0
  do i = 1,num_interfaces_crust_mantle1
    inum = 0
    do j = 1,nibool_interfaces_crust_mantle(i)
      if (isibool_interface_cm(j,i)) then
        inum = inum+1
        ibool_interfaces_crust_mantle1(inum,i) = nmir_cm(ibool_interfaces_crust_mantle(j,i))
      endif
    enddo
  enddo
  deallocate(nmir_cm)

  ! transition infinite
  if (ADD_TRINF) then
    num_interfaces_trinfinite1 = num_interfaces_trinfinite
    allocate(my_neighbors_trinfinite1(num_interfaces_trinfinite1))
    if (num_interfaces_trinfinite1 > 0) my_neighbors_trinfinite1(:) = my_neighbors_trinfinite(:)

    allocate(isibool_interface_trinf(max_nibool_interfaces_trinfinite,num_interfaces_trinfinite1), &
             nibool_interfaces_trinfinite1(num_interfaces_trinfinite1))
    isibool_interface_trinf(:,:) = .false.
    nibool_interfaces_trinfinite1(:) = 0
    do i = 1,num_interfaces_trinfinite1
      do j = 1,nibool_interfaces_trinfinite(i)
        isibool_interface_trinf(j,i)=isnode_trinf(ibool_interfaces_trinfinite(j,i))
      enddo
      nibool_interfaces_trinfinite1(i)=count(isibool_interface_trinf(:,i))
    enddo
    if (num_interfaces_trinfinite1 > 0) then
      max_nibool_interfaces_trinfinite1 = maxval(nibool_interfaces_trinfinite1)
    else
      max_nibool_interfaces_trinfinite1 = 0
    endif

    allocate(ibool_interfaces_trinfinite1(max_nibool_interfaces_trinfinite1,num_interfaces_trinfinite1))
    ibool_interfaces_trinfinite1(:,:) = 0
    do i = 1,num_interfaces_trinfinite1
      inum = 0
      do j = 1,nibool_interfaces_trinfinite(i)
        if (isibool_interface_trinf(j,i)) then
          inum = inum+1
          ibool_interfaces_trinfinite1(inum,i) = nmir_trinf(ibool_interfaces_trinfinite(j,i))
        endif
      enddo
    enddo
    deallocate(nmir_trinf)
  endif

  ! infinite
  num_interfaces_infinite1 = num_interfaces_infinite
  allocate(my_neighbors_infinite1(num_interfaces_infinite1))
  if (num_interfaces_infinite1 > 0) my_neighbors_infinite1(:) = my_neighbors_infinite(:)

  allocate(isibool_interface_inf(max_nibool_interfaces_infinite,num_interfaces_infinite1), &
           nibool_interfaces_infinite1(num_interfaces_infinite1))
  isibool_interface_inf(:,:) = .false.
  nibool_interfaces_infinite1(:) = 0
  do i = 1,num_interfaces_infinite1
    do j = 1,nibool_interfaces_infinite(i)
      isibool_interface_inf(j,i) = isnode_inf(ibool_interfaces_infinite(j,i))
    enddo
    nibool_interfaces_infinite1(i) = count(isibool_interface_inf(:,i))
  enddo
  if (num_interfaces_infinite1 > 0) then
    max_nibool_interfaces_infinite1 = maxval(nibool_interfaces_infinite1)
  else
    max_nibool_interfaces_infinite1 = 0
  endif

  allocate(ibool_interfaces_infinite1(max_nibool_interfaces_infinite1,num_interfaces_infinite1))
  ibool_interfaces_infinite1(:,:) = 0
  do i = 1,num_interfaces_infinite1
    inum = 0
    do j = 1,nibool_interfaces_infinite(i)
      if (isibool_interface_inf(j,i)) then
        inum = inum+1
        ibool_interfaces_infinite1(inum,i) = nmir_inf(ibool_interfaces_infinite(j,i))
      endif
    enddo
  enddo
  deallocate(nmir_inf)

  deallocate(isibool_interface_ic)
  deallocate(isibool_interface_oc)
  deallocate(isibool_interface_cm)
  deallocate(isibool_interface_trinf)
  deallocate(isibool_interface_inf)

  allocate(gnf1(NNDOF,nnode1),isgnf1(NNDOF,nnode1))
  gnf1 = 0
  isgnf1 = .false.

  ! activate freedoms

  ! freedoms of fictitious cube in inner core are deactivated
  do i_elmt = 1,NSPEC_INNER_CORE
    if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE) cycle

    do i = 1,NGLLCUBE_INF
      ibool = inode_elmt_ic1(i,i_elmt)
      gnf1(1,inode_ic1(ibool)) = 1
    enddo
  enddo

  ! outer core
  ! all freedoms are active
  gnf1(1,inode_oc1) = 1

  ! crust-mantle
  ! all freedoms are active
  gnf1(1,inode_cm1) = 1

  ! transition infinite
  if (ADD_TRINF) then
    ! all freedoms are active
    gnf1(1,inode_trinf1) = 1
  endif

  ! infinite element
  ! all but surface nodes are activated
  do i_elmt = 1,NSPEC_INFINITE
    do k = 1,NGLLZ_INF-1
      do j = 1,NGLLY_INF
        do i = 1,NGLLX_INF
          igll = NGLLX_INF * NGLLY_INF * (k-1) + NGLLX_INF * (j-1) + i
          ibool = inode_elmt_inf1(igll,i_elmt)
          gnf1(1,inode_inf1(ibool)) = 1
        enddo
      enddo
    enddo
  enddo

  ! copy global indices from preceding partitions
  ! inner core
  fhead = 'gdof1_ic'
  write(spm,*) i_proc
  do i = 1,num_interfaces_inner_core
    j_proc = my_neighbors_inner_core(i)
    if (j_proc < i_proc) then
      write(spn,*)j_proc
      fname='tmp_gindex3D/'//trim(fhead)//trim(adjustl(spn))//'to'//trim(adjustl(spm))
      open(10,file=fname,action='read',status='old')
      read(10,*)nibool
      allocate(gghost(NNDOF,nibool))
      read(10,*)gghost
      close(10,status='delete')

      isgnf1(:,inode_ic1(ibool_interfaces_inner_core1(1:nibool,i))) = .true.
      gnf1(:,inode_ic1(ibool_interfaces_inner_core1(1:nibool,i))) = gghost(:,:)
      deallocate(gghost)
    endif
  !    ! ndof_p2p
  !    ndof_p2p(i_proc,j_proc)=ndof_p2p(i_proc,j_proc)+NNDOF*nibool_interfaces_inner_core1(i)
  enddo

  ! outer core
  fhead = 'gdof1_oc'
  write(spm,*) i_proc
  do i = 1,num_interfaces_outer_core
    j_proc = my_neighbors_outer_core(i)
    if (j_proc < i_proc) then
      write(spn,*)j_proc
      fname='tmp_gindex3D/'//trim(fhead)//trim(adjustl(spn))//'to'//trim(adjustl(spm))
      open(10,file=fname,action='read',status='old')
      read(10,*)nibool
      allocate(gghost(NNDOF,nibool))
      read(10,*)gghost
      close(10,status='delete')

      isgnf1(:,inode_oc1(ibool_interfaces_outer_core1(1:nibool,i))) = .true.
      gnf1(:,inode_oc1(ibool_interfaces_outer_core1(1:nibool,i))) = gghost(:,:)
      deallocate(gghost)
    endif
  !    ! ndof_p2p
  !    ndof_p2p(i_proc,j_proc)=ndof_p2p(i_proc,j_proc)+NNDOF*nibool_interfaces_outer_core1(i)
  enddo

  ! crust mantle
  fhead = 'gdof1_cm'
  write(spm,*) i_proc
  do i = 1,num_interfaces_crust_mantle
    j_proc = my_neighbors_crust_mantle(i)
    if (j_proc < i_proc) then
      write(spn,*)j_proc
      fname='tmp_gindex3D/'//trim(fhead)//trim(adjustl(spn))//'to'//trim(adjustl(spm))
      open(10,file=fname,action='read',status='old')
      read(10,*)nibool
      allocate(gghost(NNDOF,nibool))
      read(10,*)gghost
      close(10,status='delete')

      isgnf1(:,inode_cm1(ibool_interfaces_crust_mantle1(1:nibool,i))) = .true.
      gnf1(:,inode_cm1(ibool_interfaces_crust_mantle1(1:nibool,i))) = gghost(:,:)
      deallocate(gghost)
    endif
  !    ! ndof_p2p
  !    ndof_p2p(i_proc,j_proc)=ndof_p2p(i_proc,j_proc)+NNDOF*nibool_interfaces_crust_mantle1(i)
  enddo

  if (ADD_TRINF) then
    ! transition infinite
    fhead = 'gdof1_trinf'
    write(spm,*)i_proc
    do i = 1,num_interfaces_trinfinite
      j_proc = my_neighbors_trinfinite(i)
      if (j_proc < i_proc) then
        write(spn,*)j_proc
        fname='tmp_gindex3D/'//trim(fhead)//trim(adjustl(spn))//'to'//trim(adjustl(spm))
        open(10,file=fname,action='read',status='old')
        read(10,*)nibool
        allocate(gghost(NNDOF,nibool))
        read(10,*)gghost
        close(10,status='delete')

        isgnf1(:,inode_trinf1(ibool_interfaces_trinfinite1(1:nibool,i))) = .true.
        gnf1(:,inode_trinf1(ibool_interfaces_trinfinite1(1:nibool,i))) = gghost(:,:)
        deallocate(gghost)
      endif
    !    ! ndof_p2p
    !    ndof_p2p(i_proc,j_proc)=ndof_p2p(i_proc,j_proc)+NNDOF*nibool_interfaces_trinfinite1(i)
    enddo
  endif

  ! infinite
  fhead = 'gdof1_inf'
  write(spm,*)i_proc
  do i = 1,num_interfaces_infinite
    j_proc = my_neighbors_infinite(i)
    if (j_proc < i_proc) then
      write(spn,*)j_proc
      fname='tmp_gindex3D/'//trim(fhead)//trim(adjustl(spn))//'to'//trim(adjustl(spm))
      open(10,file=fname,action='read',status='old')
      read(10,*)nibool
      allocate(gghost(NNDOF,nibool))
      read(10,*)gghost
      close(10,status='delete')

      isgnf1(:,inode_inf1(ibool_interfaces_infinite1(1:nibool,i))) = .true.
      gnf1(:,inode_inf1(ibool_interfaces_infinite1(1:nibool,i))) = gghost(:,:)
      deallocate(gghost)
    endif
  !    ! ndof_p2p
  !    ndof_p2p(i_proc,j_proc)=ndof_p2p(i_proc,j_proc)+NNDOF*nibool_interfaces_infinite1(i)
  enddo
  print *,'Previous largest gnf1 ID:',gnf_end1

  igdof1 = gnf_end1 ! gdof
  do i_node = 1,nnode1
    if (gnf1(1,i_node) > 0 .and. .not. isgnf1(1,i_node)) then
      isgnf1(1,i_node) = .true.
      igdof1 = igdof1+1
      gnf1(1,i_node)=igdof1
    endif
  enddo
  !neq1=inode1

  allocate(gdf_ic1(NNDOF,nnode_ic1),gdf_oc1(NNDOF,nnode_oc1), &
  gdf_cm1(NNDOF,nnode_cm1),gdf_trinf1(NNDOF,nnode_trinf1), &
  gdf_inf1(NNDOF,nnode_inf1))

  ! store gdf in a region array
  gdf_ic1(:,:) = gnf1(:,inode_ic1)
  gdf_oc1(:,:) = gnf1(:,inode_oc1)
  gdf_cm1(:,:) = gnf1(:,inode_cm1)
  gdf_trinf1(:,:) = gnf1(:,inode_trinf1)
  gdf_inf1(:,:) = gnf1(:,inode_inf1)

  deallocate(inode_ic1,inode_oc1,inode_cm1,inode_trinf1,inode_inf1)

  ! save global degrees of freedom for neighboring partitions
  ! inner core
  fhead='gdof1_ic'
  write(spm,*)i_proc
  do i = 1,num_interfaces_inner_core
    j_proc = my_neighbors_inner_core(i)
    if (j_proc > i_proc) then
      write(spn,*)j_proc
      fname='tmp_gindex3D/'//trim(fhead)//trim(adjustl(spm))//'to'//trim(adjustl(spn))
      open(10,file=fname,action='write',status='replace')
      write(10,*)nibool_interfaces_inner_core1(i)
      allocate(tmpmat(NNDOF,nibool_interfaces_inner_core1(i)))
      tmpmat = gdf_ic1(:,ibool_interfaces_inner_core1(1:nibool_interfaces_inner_core1(i),i))
      !write(10,*)gdf_ic1(:,ibool_interfaces_inner_core1(1:nibool_interfaces_inner_core1(i),i))
      write(10,*)tmpmat
      deallocate(tmpmat)
      close(10)
    endif
  enddo

  ! outer core
  fhead='gdof1_oc'
  write(spm,*)i_proc
  do i = 1,num_interfaces_outer_core
    j_proc = my_neighbors_outer_core(i)
    if (j_proc > i_proc) then
      write(spn,*)j_proc
      fname='tmp_gindex3D/'//trim(fhead)//trim(adjustl(spm))//'to'//trim(adjustl(spn))
      open(10,file=fname,action='write',status='replace')
      write(10,*)nibool_interfaces_outer_core1(i)
      allocate(tmpmat(NNDOF,nibool_interfaces_outer_core1(i)))
      tmpmat = gdf_oc1(:,ibool_interfaces_outer_core1(1:nibool_interfaces_outer_core1(i),i))
      !write(10,*)gdf_oc1(:,ibool_interfaces_outer_core1(1:nibool_interfaces_outer_core1(i),i))
      write(10,*)tmpmat
      deallocate(tmpmat)
      close(10)
    endif
  enddo

  ! crust mantle
  fhead='gdof1_cm'
  write(spm,*)i_proc
  do i = 1,num_interfaces_crust_mantle
    j_proc = my_neighbors_crust_mantle(i)
    if (j_proc > i_proc) then
      write(spn,*)j_proc
      fname='tmp_gindex3D/'//trim(fhead)//trim(adjustl(spm))//'to'//trim(adjustl(spn))
      open(10,file=fname,action='write',status='replace')
      write(10,*)nibool_interfaces_crust_mantle1(i)
      allocate(tmpmat(NNDOF,nibool_interfaces_crust_mantle1(i)))
      tmpmat = gdf_cm1(:,ibool_interfaces_crust_mantle1(1:nibool_interfaces_crust_mantle1(i),i))
      !write(10,*)gdf_cm1(:,ibool_interfaces_crust_mantle1(1:nibool_interfaces_crust_mantle1(i),i))
      write(10,*)tmpmat
      deallocate(tmpmat)
      close(10)
    endif
  enddo

  if (ADD_TRINF) then
    ! transition infinite
    fhead='gdof1_trinf'
    write(spm,*)i_proc

    do i = 1,num_interfaces_trinfinite
      j_proc = my_neighbors_trinfinite(i)
      if (j_proc > i_proc) then
        write(spn,*)j_proc
        fname='tmp_gindex3D/'//trim(fhead)//trim(adjustl(spm))//'to'//trim(adjustl(spn))
        !print *,fname
        open(10,file=fname,action='write',status='replace')
        write(10,*)nibool_interfaces_trinfinite1(i)
        allocate(tmpmat(NNDOF,nibool_interfaces_trinfinite1(i)))
        tmpmat = gdf_trinf1(:,ibool_interfaces_trinfinite1(1:nibool_interfaces_trinfinite1(i),i))
        !write(10,*)gdf_trinf1(:,ibool_interfaces_trinfinite1(1:nibool_interfaces_trinfinite1(i),i))
        write(10,*)tmpmat
        deallocate(tmpmat)
        close(10)
      endif
    enddo
  endif

  ! infinite
  fhead='gdof1_inf'
  write(spm,*) i_proc

  do i = 1,num_interfaces_infinite
    j_proc = my_neighbors_infinite(i)
    if (j_proc > i_proc) then
      write(spn,*)j_proc
      fname='tmp_gindex3D/'//trim(fhead)//trim(adjustl(spm))//'to'//trim(adjustl(spn))
      !print *,fname
      open(10,file=fname,action='write',status='replace')
      write(10,*)nibool_interfaces_infinite1(i)
      allocate(tmpmat(NNDOF,nibool_interfaces_infinite1(i)))
      tmpmat = gdf_inf1(:,ibool_interfaces_infinite1(1:nibool_interfaces_infinite1(i),i))
      !write(10,*)gdf_inf1(:,ibool_interfaces_infinite1(1:nibool_interfaces_infinite1(i),i))
      write(10,*)tmpmat
      deallocate(tmpmat)
      close(10)
    endif
  enddo

  gnf_end1 = maxval(gnf1)
  print *,'Largest gnf1 ID:',gnf_end1

  write(spm,*)i_proc

  ! file output
  ! (needed for PETSc Level-1 solver setup)
  fname='DATABASES_MPI/gdof1_proc'//trim(adjustl(spm))
  open(10,file=fname,action='write',status='replace')
  write(10,*)nnode_ic1
  write(10,*)gdf_ic1
  write(10,*)nnode_oc1
  write(10,*)gdf_oc1
  write(10,*)nnode_cm1
  write(10,*)gdf_cm1
  write(10,*)nnode_trinf1
  write(10,*)gdf_trinf1
  write(10,*)nnode_inf1
  write(10,*)gdf_inf1
  close(10)

  ! debug
  print *,'done'
  print *,'*********************************************************'

  ! deallocate variables
  deallocate(inode_ic,inode_oc,inode_cm,inode_trinf,inode_inf)
  deallocate(isnode_ic,isnode_oc,isnode_cm,isnode_trinf,isnode_inf)
  deallocate(ignode_ic,ignode_oc,ignode_cm,ignode_trinf,ignode_inf)
  deallocate(isgnode_ic,isgnode_oc,isgnode_cm,isgnode_trinf,isgnode_inf)

  ! for NGLL=5
  deallocate(gnf,isgnf)
  deallocate(gdf_ic,gdf_oc,gdf_cm,gdf_trinf,gdf_inf)

  ! deallocate arrays before re-allocating in read_mesh_databases()
  ! inner core
  deallocate(my_neighbors_inner_core,nibool_interfaces_inner_core)
  deallocate(ibool_interfaces_inner_core)
  deallocate(ibool_inner_core)
  deallocate(idoubling_inner_core)

  ! outer core
  deallocate(my_neighbors_outer_core,nibool_interfaces_outer_core)
  deallocate(ibool_interfaces_outer_core)
  deallocate(ibool_outer_core)

  ! crust/mantle
  deallocate(my_neighbors_crust_mantle,nibool_interfaces_crust_mantle)
  deallocate(ibool_interfaces_crust_mantle)
  deallocate(ibool_crust_mantle)

  ! trinfinite
  if (allocated(ibool_trinfinite)) then
    deallocate(my_neighbors_trinfinite,nibool_interfaces_trinfinite)
    deallocate(ibool_interfaces_trinfinite)
    deallocate(ibool_trinfinite)
    deallocate(ibelm_bottom_trinfinite,ibelm_top_trinfinite)
    deallocate(ibelm_xmin_trinfinite,ibelm_xmax_trinfinite)
    deallocate(ibelm_ymin_trinfinite,ibelm_ymax_trinfinite)
  endif

  ! infinite
  deallocate(my_neighbors_infinite,nibool_interfaces_infinite)
  deallocate(ibool_interfaces_infinite)
  deallocate(ibool_infinite)
  deallocate(ibelm_bottom_infinite,ibelm_top_infinite)
  deallocate(ibelm_xmin_infinite,ibelm_xmax_infinite)
  deallocate(ibelm_ymin_infinite,ibelm_ymax_infinite)

  ! coupling
  deallocate(ibelm_top_inner_core)
  deallocate(ibelm_bottom_outer_core,ibelm_top_outer_core)
  deallocate(ibelm_bottom_crust_mantle,ibelm_top_crust_mantle)

  ! for NGLL=3
  deallocate(gnf1,isgnf1)
  deallocate(gdf_ic1,gdf_oc1,gdf_cm1,gdf_trinf1,gdf_inf1)
  deallocate(my_neighbors_inner_core1,nibool_interfaces_inner_core1)
  deallocate(ibool_interfaces_inner_core1)
  !
  deallocate(my_neighbors_outer_core1,nibool_interfaces_outer_core1)
  deallocate(ibool_interfaces_outer_core1)
  !
  deallocate(my_neighbors_crust_mantle1,nibool_interfaces_crust_mantle1)
  deallocate(ibool_interfaces_crust_mantle1)
  !
  deallocate(my_neighbors_trinfinite1,nibool_interfaces_trinfinite1)
  deallocate(ibool_interfaces_trinfinite1)
  !
  deallocate(my_neighbors_infinite1,nibool_interfaces_infinite1)
  deallocate(ibool_interfaces_infinite1)

  end subroutine create_gindex_for_process
