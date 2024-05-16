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

! TODO: full gravity is not working yet, needs to fully implement solver...
#ifdef USE_PETSC_NOT_WORKING_YET

!AUTHORS:
!Hom Nath Gharti
!Stefano Zhampini
!REFERENCE:
!PETSC documentation

!-------------------------------------------------------------------------------
module solver_petsc

  !#include "petsc/finclude/petscsys.h"
  !#include "petsc/finclude/petscvec.h"
  !#include "petsc/finclude/petscvec.h90"
  !#include "petsc/finclude/petscmat.h"
  #include "petsc/finclude/petscksp.h"
  !#include "petsc/finclude/petscpc.h"

  use constants_solver, only: kreal => CUSTOM_REAL,nproc => NPROCTOT_VAL, &
    KSP_ATOL,KSP_DTOL,KSP_RTOL,KSP_MAXITER, &
    KSP_ATOL1,KSP_DTOL1,KSP_RTOL1,KSP_MAXITER1
  use specfem_par, only: myrank
  use specfem_par, only: neq,ngdof,nsparse,kmat_sparse,krow_sparse,kcol_sparse, &
    kgrow_sparse,kgcol_sparse
  use specfem_par, only: neq1,ngdof1,nsparse1,kmat_sparse1,krow_sparse1, &
    kcol_sparse1,kgrow_sparse1,kgcol_sparse1,l2gdof1

  use math_library_mpi, only: maxvec,minvec

  use specfem_par, only: num_interfaces_inner_core1, &
    max_nibool_interfaces_inner_core1, &
    my_neighbors_inner_core1,nibool_interfaces_inner_core1, &
    ibool_interfaces_inner_core1,num_interfaces_outer_core1, &
    max_nibool_interfaces_outer_core1,my_neighbors_outer_core1, &
    nibool_interfaces_outer_core1,ibool_interfaces_outer_core1, &
    num_interfaces_crust_mantle1,max_nibool_interfaces_crust_mantle1, &
    my_neighbors_crust_mantle1,nibool_interfaces_crust_mantle1, &
    ibool_interfaces_crust_mantle1,num_interfaces_trinfinite1, &
    max_nibool_interfaces_trinfinite1,my_neighbors_trinfinite1, &
    nibool_interfaces_trinfinite1,ibool_interfaces_trinfinite1, &
    num_interfaces_infinite1, &
    max_nibool_interfaces_infinite1,my_neighbors_infinite1, &
    nibool_interfaces_infinite1,ibool_interfaces_infinite1

  use petscksp

  implicit none

  PetscInt,parameter :: COMMAND = 0,CG = 1,SUPERLU = 2,MUMPS = 3 ! solver type
  PetscBool      flg,flg_ch,flg_lu,flg_ilu
  PetscInt       ival,icntl
  PetscReal      val
  ! Level-1 solver
  Vec              xvec1,bvec1,uvec1,local_vec1
  Mat              Amat1,Fmat1
  KSP              ksp1
  PC               pc1
  PetscInt         iter1
  PetscInt         solver_type1 ! solver type
  ! For communications from local to global
  VecScatter             pscat1,vscat1
  ! Stores l2g map info
  ISLocalToGlobalMapping l2gmap
  !PetscBool        flg

  ! ADJOINT SIMULATIONS
  ! FOR NOW WE ASSUME THAT THE FORWARD AND ADJOINT SIMULATIONS ARE SOLVED
  ! WITH A CG METHOD, AND THAT ONLY LEVEL 1 SOLVER IS USED
  ! Adjoint Level-1 solver
  ! notes: no b_uvec1 since it seems to just be created and destroyed
  !        no bAmat1 since we can use the same stiffness matrix
  !        I think we can probably use b_pc1 = pc1 but unsure
  Vec              b_xvec1, b_bvec1, b_local_vec1
  KSP              b_ksp1
  Mat              b_Fmat1 ! not used unless MUMPS chosen (not implemented)
  PC               b_pc1
  PetscInt         b_iter1
  PetscInt         b_solver_type1 ! solver type
  ! For communications from local to global
  VecScatter             b_pscat1, b_vscat1


  ! Level-2 solver
  Vec              xvec,bvec,uvec
  Mat              Amat
  KSP              ksp
  PC               pc
  PetscErrorCode   ierr
  PetscInt         iter
  PetscInt :: nzeros_max,nzeros_min,nzerosoff_max
  PetscInt :: ngdof_part1
  PetscInt :: ig0,ig1

contains

!===============================================================================
! Level-1 solver
!===============================================================================
  subroutine petsc_initialize1()
  use specfem_par, only: ADD_TRINF,NNDOF, SIMULATION_TYPE
  use specfem_par_innercore, only: ggdof_ic1
  use specfem_par_outercore, only: ggdof_oc1
  use specfem_par_crustmantle, only: ggdof_cm1
  use specfem_par_trinfinite, only: ggdof_trinf1
  use specfem_par_infinite, only: ggdof_inf1
  implicit none
  Vec         nnzv,nzeror_gvec1,nzeror_dvec1,nzeror_ovec1,iproc_gvec1, &
              interface_gvec1,ninterface_dvec1,ninterface_ovec1,nself_gvec1
  PetscInt :: i,istart,iend,n,n1,ncol_part1,nrow_part1
  PetscInt :: nnzmax,lsize,idxinsert(neq1),ldof1(neq1)
  PetscInt,allocatable :: nzeros(:),ig_array1(:)
  PetscScalar,allocatable :: rproc_array1(:)
  PetscScalar rval,valinsert(neq1),nnzv_v(1)
  PetscOffset nnzv_i
  PetscInt, allocatable :: nnz(:)
  IS global_is,local_is, b_global_is, b_local_is

  PetscInt :: icount,igdof,ind,maxrank0,ng,ng0,ng1,np0
  PetscInt,allocatable :: inzeror_array1(:),iproc_array1(:),nzeros_row(:)
  PetscInt,allocatable :: nnzero_diag1(:),nnzero_offdiag1(:)
  PetscInt,allocatable :: nnzero_diag1r(:),nnzero_offdiag1r(:)
  PetscScalar,pointer :: nzeror_array1(:),rproc_array(:)
  PetscScalar,pointer :: nzeror_darray1(:),nzeror_oarray1(:),rnself_array1(:)
  PetscReal :: fac_ni,max_ni,pmax,pmin,rnid,rnioffd,rnd,rnoffd,rproc,zero

  PetscInt :: ir,ic,igr,igc,ir0,ic0,igr0,igc0
  PetscInt :: nd,noffd,nid,nioffd
  PetscInt :: i_bool,i_ndof

  PetscInt :: nibool,ng_interface
  PetscInt,allocatable :: ibool_interface(:),ig_interface(:),isg_interface(:), &
                          nself_array1(:)
  PetscInt, allocatable :: ninterface_darray1(:),ninterface_oarray1(:)
  PetscScalar,allocatable :: rg_interface(:),rnself_lgarray1(:)
  PetscScalar,pointer :: rninterface_darray1(:),rninterface_oarray1(:)

  character(len=10) :: char_myrank
  character(len=60) :: outf_name


  if (myrank == 0) write(*,*)
  if (myrank == 0) write(*,*) ' ---------- Initialise PETSC: ---------- '
  if (myrank == 0) write(*,*)


  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

  ! count number of nonzeros per row
  allocate(nzeros(neq1))
  nzeros = 0
  do i = 1,nsparse1
    nzeros(krow_sparse1(i))=nzeros(krow_sparse1(i))+1
  enddo
  nzeros_max=maxvec(nzeros)
  nzeros_min=minvec(nzeros)
  nzerosoff_max = nzeros_max
  !nzeros_max=4*nzeros_max
  !nzeros=nzeros
  !nzeros=5*nzeros
  if (myrank == 0) write(*,*) 'nzeros in 1th index:',nzeros(1)
  if (myrank == 0) write(*,*) 'ngdof1:',ngdof1,' nzeros_max:',nzeros_max,' nzeros_min:', &
  nzeros_min,count(krow_sparse1 == 1)

  ! precompute ownership range OR partion layout
  !
  ng1 = ngdof1/nproc
  ng0 = ceiling(real(ngdof1)/real(nproc))

  np0 = ngdof1-nproc*ng1

  if (np0 == 0) then
  ! ng0=ng1
  ! all processors have equal gdofs
    ng = ng0
    ig0 = myrank*ng0 ! 0-based index
    ig1 = ig0+ng0-1
  else if (np0 > 0) then
  ! first np0 processors have ng0 gdofs each and remainging processors have ng1
  ! gdofs each
    maxrank0 = np0-1 ! myrank is 0-based
    if (myrank <= maxrank0) then
      ng = ng0
      ig0 = myrank*ng0 ! 0-based index
      ig1 = ig0+ng0-1
    else !myrank > maxrank0
      ng = ng1
      ig0=np0*ng0+(myrank-np0)*ng1 ! 0-based index
      ig1 = ig0+ng1-1
    endif
  else
  ! Error
    write(*,*) 'ERROR: illegal value of "np0"!'
    stop
  endif
  !if (myrank==0) write(*,*) 'OK0:',ng0,ng1,ng,ig0,ig1
  !call sync_all
  allocate(nzeros_row(ng))
  nzeros_row = 0
  !if (myrank==0) then
  !  open(1,file='test_file_proc1',action='write',status='replace')
  !  write(1,*)ng,ig0,ig1
  !endif
  do i = 1,nsparse1
    if (kgrow_sparse1(i)-1 >= ig0 .and. kgrow_sparse1(i)-1 <= ig1) then
      ind=kgrow_sparse1(i)-ig0 ! Fortran indexing
      !if (myrank==0) write(1,*) ind,kgrow_sparse1(i)
      nzeros_row(ind)=nzeros_row(ind)+1
    endif
  enddo
  !if (myrank==0) close(1)
  !nzeros_row=2*nzeros_row
  !if (myrank==0) write(*,*) 'OK1:',nzeros_row(1),minval(nzeros_row),maxval(nzeros_row)
  !call sync_all
  write(char_myrank,'(i4)')myrank
  !outf_name='precomp_nonzeros'//trim(adjustl(char_myrank))
  !open(1,file=outf_name,action='write',status='replace')
  !write(1,'(i4)')nzeros_row
  !close(1)

  call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,ngdof1,xvec1,ierr)
  CHKERRA(ierr)
  call VecDuplicate(xvec1,bvec1,ierr)
  CHKERRA(ierr)
  call VecDuplicate(xvec1,uvec1,ierr)
  CHKERRA(ierr)
  call VecDuplicate(xvec1,nzeror_gvec1,ierr)
  CHKERRA(ierr)
  call VecDuplicate(xvec1,nzeror_dvec1,ierr)
  CHKERRA(ierr)
  call VecDuplicate(xvec1,nzeror_ovec1,ierr)
  CHKERRA(ierr)
  call VecDuplicate(xvec1,iproc_gvec1,ierr)
  CHKERRA(ierr)
  call VecDuplicate(xvec1,interface_gvec1,ierr)
  CHKERRA(ierr)
  call VecDuplicate(xvec1,nself_gvec1,ierr)
  CHKERRA(ierr)
  call VecDuplicate(xvec1,ninterface_dvec1,ierr)
  CHKERRA(ierr)
  call VecDuplicate(xvec1,ninterface_ovec1,ierr)
  CHKERRA(ierr)

  if (SIMULATION_TYPE == 3) then
    ! Create backward xvector and associated b vec:
    call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,ngdof1,b_xvec1,ierr)
    CHKERRA(ierr)
    call VecDuplicate(b_xvec1,b_bvec1,ierr)
    CHKERRA(ierr)
  endif


  ! loval vector
  call VecCreateSeq(PETSC_COMM_SELF,neq1,local_vec1,ierr)
  if (SIMULATION_TYPE == 3) then
    call VecCreateSeq(PETSC_COMM_SELF,neq1,b_local_vec1,ierr)
  endif

  ! objects needed for global vector scattering to local vector
  ! create local and global IS (index set) objects from the array of local and
  ! global indices
  call ISCreateGeneral(PETSC_COMM_WORLD,neq1,l2gdof1(1:),PETSC_COPY_VALUES, &
  global_is,ierr)
  CHKERRA(ierr)
  call ISCreateStride(PETSC_COMM_SELF,neq1,0,1,local_is,ierr);
  CHKERRA(ierr)

  if (SIMULATION_TYPE == 3) then
    call ISCreateGeneral(PETSC_COMM_WORLD,neq1,l2gdof1(1:),PETSC_COPY_VALUES, &
    b_global_is,ierr)
    CHKERRA(ierr)
    call ISCreateStride(PETSC_COMM_SELF,neq1,0,1,b_local_is,ierr);
    CHKERRA(ierr)
  endif



  ! create VecScatter object which is needed to scatter PETSc parallel vectors
  call VecScatterCreate(bvec1,global_is,local_vec1,local_is,vscat1,ierr)
  CHKERRA(ierr)

  if (SIMULATION_TYPE == 3) then
    call VecScatterCreate(b_bvec1, b_global_is, b_local_vec1, b_local_is, b_vscat1, ierr)
    CHKERRA(ierr)
  endif
  call ISDestroy(global_is,ierr) ! no longer necessary
  call ISDestroy(local_is,ierr)  ! no longer necessary
  call ISDestroy(b_global_is,ierr) ! no longer necessary
  call ISDestroy(b_local_is,ierr)  ! no longer necessary


  ! assign owner processor ID to each gdof (or row)
  allocate(ig_array1(ng),rproc_array1(ng))
  ig_array1=(/ (i,i = ig0,ig1) /)
  !rproc=real(myrank)
  rproc_array1=real(myrank)
  call VecSetValues(iproc_gvec1,ng,ig_array1,rproc_array1,INSERT_VALUES,ierr);
  CHKERRA(ierr)
  deallocate(ig_array1,rproc_array1)
  call VecAssemblyBegin(iproc_gvec1,ierr)
  CHKERRA(ierr)
  call VecAssemblyEnd(iproc_gvec1,ierr)
  CHKERRA(ierr)
  call VecMin(iproc_gvec1,PETSC_NULL_INTEGER,pmin,ierr)
  call VecMax(iproc_gvec1,PETSC_NULL_INTEGER,pmax,ierr)
  if (myrank == 0) write(*,*) 'iproc range',pmin,pmax; call sync_all
  !call VecGetArrayF90(iproc_gvec1,rproc_array,ierr)
  !CHKERRA(ierr)
  !allocate(iproc_array(ng))
  !iproc_array=int(rproc_array(1:n))
  !call VecRestoreArrayF90(iproc_gvec1,rproc_array,ierr)
  !CHKERRA(ierr)
  ! copy solution to local array
  allocate(iproc_array1(neq1),rproc_array1(neq1))
  call scatter_globalvec1(iproc_gvec1, rproc_array1)
  iproc_array1=int(rproc_array1)
  if (myrank == 0) write(*,*) 'vector3 iproc',minval(iproc_array1),maxval(iproc_array1)
  call sync_all
  !!TODO: use local scatter
  !call VecScatterCreateToAll(iproc_gvec1,pscat1,iproc_garray,ierr);
  !call VecScatterBegin(pscat1,iproc_gvec1,iproc_garray,INSERT_VALUES,SCATTER_FORWARD,ierr);
  !call VecScatterEnd(pscat1,iproc_gvec1,iproc_garray,INSERT_VALUES,SCATTER_FORWARD,ierr);
  !call VecScatterDestroy(pscat1);

  ! assign interface ID to each gdofs
  rval = 1.0
  ! inner core
  do i = 1,num_interfaces_inner_core1
      nibool=nibool_interfaces_inner_core1(i)
      allocate(ibool_interface(nibool))
      ibool_interface=ibool_interfaces_inner_core1(1:nibool,i)
      !ng_interface=nibool*NNDOF
      !allocate(ig_interface(ng_interface),rg_interface(ng_interface))
      !ig_interface=reshape(ggdof_ic1(1,ibool_interface), (/ ng_interface /) )
      !ig_interface=ig_interface-1
      !rg_interface=1.0
      !call VecSetValues(interface_gvec1,ng_interface,ig_interface,rg_interface,INSERT_VALUES,ierr);
      !deallocate(ibool_interface,ig_interface,rg_interface)
      do i_bool = 1,nibool
        do i_ndof = 1,NNDOF
          igdof = ggdof_ic1(i_ndof,ibool_interface(i_bool))-1
          if (igdof >= 0) call VecSetValues(interface_gvec1,1,igdof,rval, &
          INSERT_VALUES,ierr);
        enddo
      enddo
      deallocate(ibool_interface)
  enddo
  !if (myrank==0) write(*,*) 'OK'
  !call sync_all
  !! stop all the MPI processes, and exit
  !call MPI_FINALIZE(ierr)

  ! outer core
  do i = 1,num_interfaces_outer_core1
      nibool=nibool_interfaces_outer_core1(i)
      allocate(ibool_interface(nibool))
      ibool_interface=ibool_interfaces_outer_core1(1:nibool,i)
      !ng_interface=nibool*NNDOF
      !allocate(ig_interface(ng_interface),rg_interface(ng_interface))
      !ig_interface=reshape(ggdof_ic1(1,ibool_interface), (/ ng_interface /) )
      !ig_interface=ig_interface-1
      !rg_interface=1.0
      !call VecSetValues(interface_gvec1,ng_interface,ig_interface,rg_interface,INSERT_VALUES,ierr);
      !deallocate(ibool_interface,ig_interface,rg_interface)
      do i_bool = 1,nibool
        do i_ndof = 1,NNDOF
          igdof = ggdof_oc1(i_ndof,ibool_interface(i_bool))-1
          if (igdof >= 0) call VecSetValues(interface_gvec1,1,igdof,rval, &
          INSERT_VALUES,ierr);
        enddo
      enddo
      deallocate(ibool_interface)
  enddo

  ! crust mantle
  do i = 1,num_interfaces_crust_mantle1
      nibool=nibool_interfaces_crust_mantle1(i)
      allocate(ibool_interface(nibool))
      ibool_interface=ibool_interfaces_crust_mantle1(1:nibool,i)
      !ng_interface=nibool*NNDOF
      !allocate(ig_interface(ng_interface),rg_interface(ng_interface))
      !ig_interface=reshape(ggdof_ic1(1,ibool_interface), (/ ng_interface /) )
      !ig_interface=ig_interface-1
      !rg_interface=1.0
      !call VecSetValues(interface_gvec1,ng_interface,ig_interface,rg_interface,INSERT_VALUES,ierr);
      !deallocate(ibool_interface,ig_interface,rg_interface)
      do i_bool = 1,nibool
        do i_ndof = 1,NNDOF
          igdof = ggdof_cm1(i_ndof,ibool_interface(i_bool))-1
          if (igdof >= 0) call VecSetValues(interface_gvec1,1,igdof,rval, &
          INSERT_VALUES,ierr);
        enddo
      enddo
      deallocate(ibool_interface)
  enddo

  ! transition infinite
  if (ADD_TRINF) then
  do i = 1,num_interfaces_trinfinite1
      nibool=nibool_interfaces_trinfinite1(i)
      allocate(ibool_interface(nibool))
      ibool_interface=ibool_interfaces_trinfinite1(1:nibool,i)
      !ng_interface=nibool*NNDOF
      !allocate(ig_interface(ng_interface),rg_interface(ng_interface))
      !ig_interface=reshape(ggdof_ic1(1,ibool_interface), (/ ng_interface /) )
      !ig_interface=ig_interface-1
      !rg_interface=1.0
      !call VecSetValues(interface_gvec1,ng_interface,ig_interface,rg_interface,INSERT_VALUES,ierr);
      !deallocate(ibool_interface,ig_interface,rg_interface)
      do i_bool = 1,nibool
        do i_ndof = 1,NNDOF
          igdof = ggdof_trinf1(i_ndof,ibool_interface(i_bool))-1
          if (igdof >= 0) call VecSetValues(interface_gvec1,1,igdof,rval, &
          INSERT_VALUES,ierr);
        enddo
      enddo
      deallocate(ibool_interface)
  enddo
  endif

  ! infinite
  do i = 1,num_interfaces_infinite1
      nibool=nibool_interfaces_infinite1(i)
      allocate(ibool_interface(nibool))
      ibool_interface=ibool_interfaces_infinite1(1:nibool,i)
      !ng_interface=nibool*NNDOF
      !allocate(ig_interface(ng_interface),rg_interface(ng_interface))
      !ig_interface=reshape(ggdof_ic1(1,ibool_interface), (/ ng_interface /) )
      !ig_interface=ig_interface-1
      !rg_interface=1.0
      !call VecSetValues(interface_gvec1,ng_interface,ig_interface,rg_interface,INSERT_VALUES,ierr);
      !deallocate(ibool_interface,ig_interface,rg_interface)
      do i_bool = 1,nibool
        do i_ndof = 1,NNDOF
          igdof = ggdof_inf1(i_ndof,ibool_interface(i_bool))-1
          if (igdof >= 0) call VecSetValues(interface_gvec1,1,igdof,rval, &
          INSERT_VALUES,ierr);
        enddo
      enddo
      deallocate(ibool_interface)
  enddo

  call VecAssemblyBegin(interface_gvec1,ierr)
  CHKERRA(ierr)
  call VecAssemblyEnd(interface_gvec1,ierr)
  CHKERRA(ierr)

  !call sync_all
  !! stop all the MPI processes, and exit
  !call MPI_FINALIZE(ierr)

  ! copy solution to local array
  allocate(isg_interface(neq1),rg_interface(neq1))
  call scatter_globalvec1(interface_gvec1,rg_interface)
  isg_interface=int(rg_interface)

  ! estimate correction for the number of nonzero entries in the diagonal and
  ! nondiagonal portion
  ! self interface
  !rval=-1.0
  !call VecSet(nself_gvec1,rval,ierr) ! subtract self
  rval = 1.0
  do i = 1,neq1
    if (isg_interface(i) == 1) then
      call VecSetValues(nself_gvec1,1,l2gdof1(i),rval,ADD_VALUES,ierr);
    endif
  enddo
  call VecAssemblyBegin(nself_gvec1,ierr)
  CHKERRA(ierr)
  call VecAssemblyEnd(nself_gvec1,ierr)
  CHKERRA(ierr)
  call VecGetLocalSize(nself_gvec1,n,ierr)

  allocate(rnself_lgarray1(neq1))
  call scatter_globalvec1(nself_gvec1, rnself_lgarray1)
  call VecGetArrayF90(nself_gvec1,rnself_array1,ierr)
  allocate(nself_array1(n))
  nself_array1 = int(rnself_array1(1:n))
  where(nself_array1 > 0)nself_array1=nself_array1-1 ! subtract self
  call VecRestoreArrayF90(nself_gvec1,rnself_array1,ierr)
  call VecDestroy(nself_gvec1,ierr)

  if (myrank == 0) write(*,*) 'maximum value of nself:',maxval(nself_array1)
  !call sync_all
  !! stop all the MPI processes, and exit
  !call MPI_FINALIZE(ierr)
  !stop
  !! count nonzero entries in the diagonal and nondiagonal portion
  !zero=0.
  !call VecSet(nzeror_dvec1,zero,ierr)
  !call VecSet(nzeror_ovec1,zero,ierr)

  !outf_name='isg_interface'//trim(adjustl(char_myrank))
  !open(1,file=outf_name,action='write',status='replace')
  !write(1,'(i4)')isg_interface
  !close(1)

  !if (myrank==0) open(11,file='test_interface',action='write',status='replace')
  ! factor for maximum number of interfaces for each nondiagonal entry of the
  ! stiffness matrix
  ! the factor below is valid ONLY for rectagular partitioning of the global model
  max_ni = 8.0
  fac_ni = 0.0

  ! first element
  igr0=kgrow_sparse1(1)-1
  igc0=kgcol_sparse1(1)-1
  ir0=krow_sparse1(1)
  ic0=kcol_sparse1(1)
  nd = 0; noffd = 0
  rnid = 0.; rnioffd = 0.
  if (iproc_array1(ir0) == iproc_array1(ic0)) then
    nd = 1;
    if (igr0 /= igc0 .and. rnself_lgarray1(ir0) > 1.0 .and. rnself_lgarray1(ic0) > 1.0) then
      fac_ni = min(max_ni,min(rnself_lgarray1(ir0),rnself_lgarray1(ic0)))
      rnid = 1.0/fac_ni
    endif
  else
    noffd = 1
    if (igr0 /= igc0 .and. rnself_lgarray1(ir0) > 1.0 .and. rnself_lgarray1(ic0) > 1.0) then
      fac_ni = min(max_ni,min(rnself_lgarray1(ir0),rnself_lgarray1(ic0)))
      rnioffd = 1.0/fac_ni
    endif
  endif
  do i = 2,nsparse1
    igr=kgrow_sparse1(i)-1
    igc=kgcol_sparse1(i)-1
    ir=krow_sparse1(i)
    ic=kcol_sparse1(i)
    if (l2gdof1(ir) /= igr.or.l2gdof1(ic) /= igc) then
      write(*,*) 'strange:',l2gdof1(ir),igr,l2gdof1(ic),igc
      stop
    endif
    if (igr /= igr0) then
      ! new row starts
      ! set values computed so far
      rnd=real(nd)
      rnoffd=real(noffd)
      call VecSetValues(nzeror_dvec1,1,igr0,rnd,ADD_VALUES,ierr)
      CHKERRA(ierr)
      call VecSetValues(nzeror_ovec1,1,igr0,rnoffd,ADD_VALUES,ierr)
      CHKERRA(ierr)

      call VecSetValues(ninterface_dvec1,1,igr0,rnid,ADD_VALUES,ierr)
      CHKERRA(ierr)
      call VecSetValues(ninterface_ovec1,1,igr0,rnioffd,ADD_VALUES,ierr)
      CHKERRA(ierr)

      ! reset
      nd = 0; noffd = 0
      rnid = 0.; rnioffd = 0.
      igr0=igr !kgrow_sparse1(i)-1
      igc0=igc !kgcol_sparse1(i)-1
      ir0=ir !krow_sparse1(i)
      ic0=ic !kcol_sparse1(i)

      if (iproc_array1(ir0) == iproc_array1(ic0)) then
        nd = 1;
        if (igr0 /= igc0 .and. rnself_lgarray1(ir0) > 0.0 .and. rnself_lgarray1(ic0) > 0.0) then
          fac_ni = min(max_ni,min(rnself_lgarray1(ir0),rnself_lgarray1(ic0)))
          rnid = 1.0/fac_ni
        endif
      else
        noffd = 1
        if (igr0 /= igc0 .and. rnself_lgarray1(ir0) > 0.0 .and. rnself_lgarray1(ic0) > 0.0) then
          fac_ni = min(max_ni,min(rnself_lgarray1(ir0),rnself_lgarray1(ic0)))
          rnioffd = 1.0/fac_ni
        endif
      endif
    else
      !if (myrank==0) write(11,*) ir,ic,isg_interface(ir),isg_interface(ic)
      ! count
      if (iproc_array1(ir) == iproc_array1(ic)) then
        nd = nd+1;
        if (igr /= igc .and. rnself_lgarray1(ir) > 0.0 .and. rnself_lgarray1(ic) > 0.0) then
          fac_ni = min(max_ni,min(rnself_lgarray1(ir),rnself_lgarray1(ic)))
          rnid=rnid+(1.0/fac_ni)
        endif
      else
        noffd = noffd+1
        if (igr /= igc .and. rnself_lgarray1(ir) > 0.0 .and. rnself_lgarray1(ic) > 0.0) then
          fac_ni = min(max_ni,min(rnself_lgarray1(ir),rnself_lgarray1(ic)))
          rnioffd=rnioffd+(1.0/fac_ni)
        endif
      endif
    endif
    if (i == nsparse1) then
      ! for last
      rnd=real(nd)
      rnoffd=real(noffd)
      call VecSetValues(nzeror_dvec1,1,igr0,rnd,ADD_VALUES,ierr)
      CHKERRA(ierr)
      call VecSetValues(nzeror_ovec1,1,igr0,rnoffd,ADD_VALUES,ierr)
      CHKERRA(ierr)

      call VecSetValues(ninterface_dvec1,1,igr0,rnid,ADD_VALUES,ierr)
      CHKERRA(ierr)
      call VecSetValues(ninterface_ovec1,1,igr0,rnioffd,ADD_VALUES,ierr)
      CHKERRA(ierr)
    endif
  enddo
  deallocate(krow_sparse1,kcol_sparse1)
  ! Assemble vectors globally
  call VecAssemblyBegin(nzeror_dvec1,ierr)
  call VecAssemblyEnd(nzeror_dvec1,ierr)
  call VecAssemblyBegin(nzeror_ovec1,ierr)
  call VecAssemblyEnd(nzeror_ovec1,ierr)

  call VecAssemblyBegin(ninterface_dvec1,ierr)
  call VecAssemblyEnd(ninterface_dvec1,ierr)
  call VecAssemblyBegin(ninterface_ovec1,ierr)
  call VecAssemblyEnd(ninterface_ovec1,ierr)

  ! apply correction for repeatition due to interfaces
  call VecGetLocalSize(nzeror_dvec1,n,ierr)
  call VecGetArrayF90(nzeror_dvec1,nzeror_darray1,ierr)
  allocate(nnzero_diag1(n))
  nnzero_diag1 = int(nzeror_darray1(1:n))
  nnzero_diag1 = nnzero_diag1-nself_array1

  if (myrank == 0) write(*,*) n,minval(nzeror_darray1),maxval(nzeror_darray1), &
  minval(nnzero_diag1),maxval(nnzero_diag1)
  call sync_all

  call VecRestoreArrayF90(nzeror_dvec1,nzeror_darray1,ierr)
  call VecDestroy(nzeror_dvec1,ierr)

  call VecGetArrayF90(nzeror_ovec1,nzeror_oarray1,ierr)
  allocate(nnzero_offdiag1(n))
  nnzero_offdiag1 = int(nzeror_oarray1(1:n))

  call VecRestoreArrayF90(nzeror_ovec1,nzeror_oarray1,ierr)
  call VecDestroy(nzeror_ovec1,ierr)

  ! correction
  ! I do not know why but there are some DOFs where the correction exceeds by 4 or
  ! 8 therefore to be safe we need to subtract this from all
  call VecGetArrayF90(ninterface_dvec1,rninterface_darray1,ierr)
  !where(rninterface_darray1>0.0 .and. rninterface_darray1 < 1.0)rninterface_darray1=1.0
  allocate(ninterface_darray1(n))
  ninterface_darray1 = int(rninterface_darray1(1:n))
  call VecRestoreArrayF90(ninterface_dvec1,rninterface_darray1,ierr)
  call VecDestroy(ninterface_dvec1,ierr)
  where(ninterface_darray1 > 0)ninterface_darray1=ninterface_darray1-4
  where(ninterface_darray1 < 0)ninterface_darray1=0


  call VecGetArrayF90(ninterface_ovec1,rninterface_oarray1,ierr)
  !where(rninterface_oarray1>0.0 .and. rninterface_oarray1 < 1.0)rninterface_oarray1=1.0
  allocate(ninterface_oarray1(n))
  ninterface_oarray1 = int(rninterface_oarray1(1:n))
  call VecRestoreArrayF90(ninterface_ovec1,rninterface_oarray1,ierr)
  call VecDestroy(ninterface_ovec1,ierr)
  where(ninterface_oarray1 > 0)ninterface_oarray1=ninterface_oarray1-8
  where(ninterface_oarray1 < 0)ninterface_oarray1=0


  nnzero_diag1 = nnzero_diag1-ninterface_darray1
  nnzero_offdiag1 = nnzero_offdiag1-ninterface_oarray1

  do i = 1,nsparse1
    rval = 1.
    igdof=kgrow_sparse1(i)-1 ! Fortran index
    call VecSetValues(nzeror_gvec1,1,igdof,rval,ADD_VALUES,ierr);
    CHKERRA(ierr)
  enddo
  call VecAssemblyBegin(nzeror_gvec1,ierr)
  CHKERRA(ierr)
  call VecAssemblyEnd(nzeror_gvec1,ierr)
  CHKERRA(ierr)
  call VecGetLocalSize(nzeror_gvec1,n,ierr)
  CHKERRA(ierr)
  if (myrank == 0) write(*,*) 'size of vector:',ng,n,minval(kgrow_sparse1),ig0
  call VecGetArrayF90(nzeror_gvec1,nzeror_array1,ierr)
  CHKERRA(ierr)
  allocate(inzeror_array1(n))
  inzeror_array1 = int(nzeror_array1(1:n))
  call VecRestoreArrayF90(nzeror_gvec1,nzeror_array1,ierr)
  CHKERRA(ierr)
  call VecDestroy(nzeror_gvec1,ierr)
  CHKERRA(ierr)


  ! Create the stiffness matrix (same for forward/adjoint simulations)
  call MatCreate(PETSC_COMM_WORLD,Amat1,ierr)
  call MatSetType(Amat1,MATMPIAIJ,ierr)
  CHKERRA(ierr)
  call MatSetSizes(Amat1,PETSC_DECIDE,PETSC_DECIDE,ngdof1,ngdof1,ierr)
  CHKERRA(ierr)

  call MatMPIAIJSetPreallocation(Amat1,nzeros_max,inzeror_array1,nzeros_max, &
  inzeror_array1,ierr)
  CHKERRA(ierr)

  call MatSetFromOptions(Amat1,ierr)
  CHKERRA(ierr)

  call MatGetOwnershipRange(Amat1,istart,iend,ierr)
  CHKERRA(ierr)
  call sync_all
  if (istart /= ig0 .or. iend-1 /= ig1) then
    write(*,*) 'ERROR: ownership range mismatch!'
    write(*,*) 'ownership range:',myrank,istart,ig0,iend-1,ig1,nzeros_row(1)
    stop
  endif
  deallocate(nzeros)

  ! Create forward solver
  !solver_type1=CG
  !call create_linear_solver(solver_type1, ksp1, Amat1, pc1, Fmat1)

  !call KSPSetTolerances(ksp1,KSP_RTOL1,KSP_ATOL1,KSP_DTOL1,KSP_MAXITER1,ierr)
  !CHKERRA(ierr)

  !  Set runtime options, e.g.,
  !    -ksp_type < type> -pc_type < type> -ksp_monitor -ksp_KSP_RTOL < KSP_RTOL>
  !  These options will override those specified above as long as
  !  KSPSetFromOptions() is called _after_ any other customization
  !  routines.
  !call KSPSetFromOptions(ksp1,ierr)


  ! Create adjoint solver:
  !if (SIMULATION_TYPE==3) then
  !  call create_linear_solver(solver_type1, b_ksp1, Amat1, pc1, Fmat1)
  !  call KSPSetTolerances(b_ksp1,KSP_RTOL1,KSP_ATOL1,KSP_DTOL1,KSP_MAXITER1,ierr)
  !  CHKERRA(ierr)
  !  call KSPSetFromOptions(b_ksp1,ierr)
  !  if (myrank==0) then
  !    write(*,*) ' Created adjoint linear KSP solver...'
  !  endif
  !endif




  !-------------------------------------------------------------------------------
  ! Create the linear solver and set various options
  !-------------------------------------------------------------------------------
  ! define solver type
  ! COMMAND: define from the command
  ! SUPERLU: SuperLU solver
  ! MUMPS: MUMPS solver
  !solver_type1=CG
  ! Create linear solver context

  call KSPCreate(PETSC_COMM_WORLD,ksp1,ierr)
  call KSPSetOperators(ksp1,Amat1,Amat1,ierr) ! version >= 3.5

  if (SIMULATION_TYPE == 3) then
    call KSPSetInitialGuessNonzero(ksp1,PETSC_FALSE,ierr)
  else
    call KSPSetInitialGuessNonzero(ksp1,PETSC_TRUE,ierr)
  endif

  CHKERRA(ierr)
  call KSPSetDiagonalScale(ksp1,PETSC_TRUE,ierr)
  CHKERRA(ierr)
  call KSPSetReusePreconditioner(ksp1,PETSC_TRUE,ierr)
  call KSPSetType(ksp1,KSPCG,ierr);
  CHKERRA(ierr)
  call KSPGetPC(ksp1,pc1,ierr)
  CHKERRA(ierr)
  call PCFactorSetShiftType(pc1,MAT_SHIFT_POSITIVE_DEFINITE,ierr)
  CHKERRA(ierr)
  call KSPSetTolerances(ksp1,KSP_RTOL1,KSP_ATOL1,KSP_DTOL1,KSP_MAXITER1,ierr)
  CHKERRA(ierr)
  call KSPSetFromOptions(ksp1,ierr)

  ! BACKWARD SOLVER
  call KSPCreate(PETSC_COMM_WORLD,b_ksp1,ierr)
  call KSPSetOperators(b_ksp1,Amat1,Amat1,ierr) ! version >= 3.5
  call KSPSetInitialGuessNonzero(b_ksp1,PETSC_FALSE,ierr)
  CHKERRA(ierr)
  call KSPSetDiagonalScale(b_ksp1,PETSC_TRUE,ierr)
  CHKERRA(ierr)
  call KSPSetReusePreconditioner(b_ksp1,PETSC_TRUE,ierr)
  call KSPSetType(b_ksp1,KSPCG,ierr);
  CHKERRA(ierr)
  call KSPGetPC(b_ksp1,b_pc1,ierr)
  CHKERRA(ierr)
  call PCFactorSetShiftType(b_pc1,MAT_SHIFT_POSITIVE_DEFINITE,ierr)
  CHKERRA(ierr)
  call KSPSetTolerances(b_ksp1,KSP_RTOL1,KSP_ATOL1,KSP_DTOL1,KSP_MAXITER1,ierr)
  CHKERRA(ierr)
  call KSPSetFromOptions(b_ksp1,ierr)

  if (myrank == 0) then
    write(*,*) ' ---------- Finished PETSC initialisation ---------- '
  endif

  end subroutine petsc_initialize1

!
!===============================================================================
!

  subroutine create_linear_solver(stype, l_ksp, l_Amat, l_pc, l_fmat)

  ! Create the linear solver and set various options
  ! stype: Solver type - options available are
  !   COMMAND   define from the command
  !   SUPERLU   SuperLU solver
  !   MUMPS     MUMPS solver
  ! l_ksp: the local KSP (ksp1 or b_ksp1 etc)
  ! l_Amat: local A matrix - I think always Amat1
  ! l_pc:   local preconditioner e.g pc1
  use specfem_par, only: myrank, SIMULATION_TYPE
  implicit none



  PetscInt     stype
  KSP          l_ksp
  Mat          l_Amat, l_fmat
  PC           l_pc


  ! Create linear solver context
  call KSPCreate(PETSC_COMM_WORLD,l_ksp,ierr)
  ! Set operators. Here the matrix that defines the linear system
  ! also serves as the preconditioning matrix.
  !call KSPSetOperators(ksp1,Amat1,Amat1,SAME_PRECONDITIONER,ierr) ! version < 3.5
  call KSPSetOperators(l_ksp,l_Amat,l_Amat,ierr) ! version >= 3.5

  call KSPSetInitialGuessNonzero(l_ksp,PETSC_TRUE,ierr)
  CHKERRA(ierr)
  !since the euqutions are nondimensionalized, the scaling is unnecessary?
  call KSPSetDiagonalScale(l_ksp,PETSC_TRUE,ierr)
  CHKERRA(ierr)
  call KSPSetReusePreconditioner(l_ksp,PETSC_TRUE,ierr)

  if (stype == COMMAND) then
    if (myrank == 0) write(*,*) 'Solver type: provided via command'
  else if (stype == CG) then
    ! CONJUGATE GRADIENT
    if (myrank == 0) write(*,*) 'Solver type: CG'
    call KSPSetType(l_ksp,KSPCG,ierr);
    CHKERRA(ierr)
    ! Fetch preconditioner
    call KSPGetPC(l_ksp,l_pc,ierr)
    CHKERRA(ierr)
    call PCFactorSetShiftType(l_pc,MAT_SHIFT_POSITIVE_DEFINITE,ierr)
    CHKERRA(ierr)
  else if (stype == SUPERLU) then
    ! SUPER LU
    if (myrank == 0) write(*,*) 'Solver type: SUPERLU'
    if (SIMULATION_TYPE == 3) then
      write(*,*) ' ERROR: SUPERLU not implemented for adjoint sims yet.'
      stop
    endif
    flg_ilu    = PETSC_FALSE;
    flg_lu     = PETSC_FALSE;
    ! version < 3.8.0
      ! call
      ! PetscOptionsGetBool(PETSC_NULL_CHARACTER,"-use_superlu_lu",flg_lu,flg,ierr);
        call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
          "-use_superlu_lu",flg_lu,flg,ierr);
    CHKERRA(ierr)
    !PetscOptionsGetBool(PETSC_NULL_CHARACTER,"-use_superlu_ilu",flg_ilu,flg,ierr);
    if (flg_lu .or. flg_ilu) then
      call KSPSetType(l_ksp,KSPPREONLY,ierr);
      CHKERRA(ierr)
      call KSPGetPC(l_ksp,l_pc,ierr);
      CHKERRA(ierr)
      if (flg_lu) then
        call PCSetType(l_pc,PCLU,ierr);
        CHKERRA(ierr)
      else if (flg_ilu) then
        call PCSetType(l_pc,PCILU,ierr);
        CHKERRA(ierr)
      endif
      call PCFactorSetShiftType(l_pc,MAT_SHIFT_POSITIVE_DEFINITE,ierr)
      CHKERRA(ierr)
      ! version < 3.9
      !call PCFactorSetMatSolverPackage(l_pc,MATSOLVERSUPERLU,ierr);
      call PCFactorSetMatSolverType(l_pc,MATSOLVERSUPERLU,ierr);
      CHKERRA(ierr)
      ! version < 3.9
      !call PCFactorSetUpMatSolverPackage(l_pc,ierr); ! call MatGetFactor() to create F
      call PCFactorSetUpMatSolverType(l_pc,ierr); ! call MatGetFactor() to create F
      CHKERRA(ierr)
      call PCFactorGetMatrix(l_pc,l_fmat,ierr);
      CHKERRA(ierr)
      !call MatSuperluSetILUDropTol(l_fmat,1.e-8,ierr);
      !CHKERRA(ierr)
    endif
  else if (stype == MUMPS) then
    if (myrank == 0) write(*,*) 'Solver type: MUMPS'
    write(*,*) 'ERROR - WE commented out MUMPS stuff due to syntax error'
    stop
    flg_lu    = PETSC_FALSE;
    flg_ch    = PETSC_FALSE;
    ! version < 3.8.0
    ! call PetscOptionsGetBool(PETSC_NULL_CHARACTER,"-use_mumps_ch",flg_ch,flg,ierr);
     call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
          "-use_mumps_ch",flg_ch,flg,ierr);
    if (flg_lu .or. flg_ch) then
      call KSPSetType(l_ksp,KSPPREONLY,ierr);
      call KSPGetPC(l_ksp,l_pc,ierr);
      if (flg_lu) then
        call PCSetType(l_pc,PCLU,ierr);
      else if (flg_ch) then
        call MatSetOption(l_Amat,MAT_SPD,PETSC_TRUE,ierr); ! set MUMPS id%SYM=1
        call PCSetType(l_pc,PCCHOLESKY,ierr);
      endif
      call PCFactorSetShiftType(l_pc,MAT_SHIFT_POSITIVE_DEFINITE,ierr)
      CHKERRA(ierr)
      ! version < 3.9
      !call PCFactorSetMatSolverPackage(l_pc,MATSOLVERMUMPS,ierr);
      !call PCFactorSetUpMatSolverPackage(l_pc,ierr); ! call MatGetFactor() to create F
      !call PCFactorSetMatSolverType(l_pc,MATSOLVERMUMPS,ierr);
      !call PCFactorSetUpMatSolverType(l_pc,ierr); ! call MatGetFactor() to create F
      !call PCFactorGetMatrix(l_pc,l_fmat,ierr);
      icntl = 7; ival = 2;
      !call MatMumpsSetIcntl(l_fmat,icntl,ival,ierr);
      icntl = 1; val = 0.0;
      !call MatMumpsSetCntl(l_fmat,icntl,val,ierr);
    endif
  endif

  end subroutine create_linear_solver

!
!===============================================================================
!

  subroutine petsc_set_matrix1()

  use math_library_mpi, only: maxscal,minscal
  use specfem_par, only: IFLAG_IN_FICTITIOUS_CUBE,NSPEC_INNER_CORE, &
  NSPEC_OUTER_CORE,NSPEC_CRUST_MANTLE,NSPEC_TRINFINITE,NSPEC_INFINITE
  use specfem_par, only: NEDOF
  use specfem_par, only: NGLLCUBE_INF,NEDOF1
  use specfem_par_innercore, only: ggdof_ic1,storekmat_inner_core1, &
  idoubling_inner_core,inode_elmt_ic1
  use specfem_par_outercore, only: ggdof_oc1,storekmat_outer_core1,inode_elmt_oc1
  use specfem_par_crustmantle, only: ggdof_cm1,storekmat_crust_mantle1, &
  inode_elmt_cm1
  use specfem_par_trinfinite, only: ggdof_trinf1,storekmat_trinfinite1, &
  inode_elmt_trinf1
  use specfem_par_infinite, only: ggdof_inf1,storekmat_infinite1,inode_elmt_inf1

  implicit none
  integer :: i,i_elmt,j,ncount
  integer :: ggdof_elmt(NEDOF1),idof(NEDOF1),igdof(NEDOF1)

  PetscInt irow,istart,iend,ndiag,noffdiag
  integer :: ncols,ncols_val
  integer,allocatable :: cols(:)
  real(kind=8),allocatable :: vals(:)

  character(len=10) :: char_myrank
  character(len=60) :: outf_name
  Vec Adiag1,lAdiag1
  PetscScalar,pointer :: arrayAdiag1(:)
  PetscReal :: maxdiag1,mindiag1
  PetscInt n

  call MatZeroEntries(Amat1,ierr)
  CHKERRA(ierr)

  ! inner core
  do i_elmt = 1,NSPEC_INNER_CORE
    if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE)cycle
    ggdof_elmt = reshape(ggdof_ic1(:,inode_elmt_ic1(:,i_elmt)),(/NEDOF1/))
    ggdof_elmt = ggdof_elmt-1 ! petsc index starts from 0
    ncount = 0; idof=-1; igdof=-1
    do i = 1,NEDOF1
      do j = 1,NEDOF1
      if (ggdof_elmt(i) >= 0.and.ggdof_elmt(j) >= 0.and.                          &
      storekmat_inner_core1(i,j,i_elmt) /= 0.0_kreal) then
        !ncount=ncount+1
        !idof(ncount)=i
        !igdof(ncount)=ggdof_elmt(i)
        call MatSetValues(Amat1,1,ggdof_elmt(i),1,ggdof_elmt(j), &
        storekmat_inner_core1(i,j,i_elmt),ADD_VALUES,ierr)
        CHKERRA(ierr)
      endif
      enddo
    enddo
    !call MatSetValues(Amat,ncount,igdof(1:ncount),ncount,igdof(1:ncount), &
    !storekmat_inner_core(idof(1:ncount),idof(1:ncount),i_elmt),ADD_VALUES,ierr)
    !CHKERRA(ierr)
  enddo
  !! inner core
  !do i_elmt=1,NSPEC_INNER_CORE
  !  if (idoubling_inner_core(i_elmt)==IFLAG_IN_FICTITIOUS_CUBE)cycle
  !  ggdof_elmt=reshape(ggdof_ic(:,inode_elmt_ic(:,i_elmt)),(/NEDOF/))
  !  ggdof_elmt=ggdof_elmt-1 ! petsc index starts from 0
  !  ncount=0; idof=-1; igdof=-1
  !  do i=1,NEDOF
  !    if (ggdof_elmt(i) >= 0) then
  !      ncount=ncount+1
  !      idof(ncount)=i
  !      igdof(ncount)=ggdof_elmt(i)
  !    endif
  !  enddo
  !  call MatSetValues(Amat,ncount,igdof(1:ncount),ncount,igdof(1:ncount), &
  !  storekmat_inner_core(idof(1:ncount),idof(1:ncount),i_elmt),ADD_VALUES,ierr)
  !  CHKERRA(ierr)
  !enddo
  deallocate(storekmat_inner_core1)
  !if (myrank==0) write(*,*) 'IC kmat done1!'; call sync_all
  !if (myrank==0) then
  !  open(1111,file='debug.log',action='write')
  !  write(1111,*)ggdof_oc1
  !  !do i_elmt=1,NSPEC_OUTER_CORE
  !  !  ggdof_elmt=reshape(ggdof_oc1(:,inode_elmt_oc1(:,i_elmt)),(/NEDOF1/))
  !  !  ggdof_elmt=ggdof_elmt-1 ! petsc index starts from 0
  !  !  ncount=0; idof=-1; igdof=-1
  !  !  do i=1,NEDOF1
  !  !    do j=1,NEDOF1
  !  !    if (ggdof_elmt(i) >= 0.and.ggdof_elmt(j) >= 0.and.                          &
  !  !    storekmat_outer_core1(i,j,i_elmt) /= 0.0_kreal) then
  !  !      write(1111,*)ggdof_elmt(i),ggdof_elmt(j)
  !  !    endif
  !  !    enddo
  !  !  enddo
  !  !enddo
  !  close(1111)
  !endif
  !call sync_all

  ! outer core
  do i_elmt = 1,NSPEC_OUTER_CORE
    ggdof_elmt = reshape(ggdof_oc1(:,inode_elmt_oc1(:,i_elmt)),(/NEDOF1/))
    ggdof_elmt = ggdof_elmt-1 ! petsc index starts from 0
    ncount = 0; idof=-1; igdof=-1
    do i = 1,NEDOF1
      do j = 1,NEDOF1
      if (ggdof_elmt(i) >= 0.and.ggdof_elmt(j) >= 0.and.                          &
      storekmat_outer_core1(i,j,i_elmt) /= 0.0_kreal) then
        call MatSetValues(Amat1,1,ggdof_elmt(i),1,ggdof_elmt(j), &
        storekmat_outer_core1(i,j,i_elmt),ADD_VALUES,ierr)
        CHKERRA(ierr)
      endif
      enddo
    enddo
  enddo
  !! outer core
  !do i_elmt=1,NSPEC_OUTER_CORE
  !  ggdof_elmt=reshape(ggdof_oc(:,inode_elmt_oc(:,i_elmt)),(/NEDOF/))
  !  ggdof_elmt=ggdof_elmt-1 ! petsc index starts from 0
  !  ncount=0; idof=-1; igdof=-1
  !  do i=1,NEDOF
  !    if (ggdof_elmt(i) >= 0) then
  !      ncount=ncount+1
  !      idof(ncount)=i
  !      igdof(ncount)=ggdof_elmt(i)
  !    endif
  !  enddo
  !  call MatSetValues(Amat,ncount,igdof(1:ncount),ncount,igdof(1:ncount), &
  !  storekmat_outer_core(idof(1:ncount),idof(1:ncount),i_elmt),ADD_VALUES,ierr)
  !  CHKERRA(ierr)
  !enddo
  deallocate(storekmat_outer_core1)
  !if (myrank==0) write(*,*) 'OC kmat done1!'; call sync_all
  ! crust mantle
  do i_elmt = 1,NSPEC_CRUST_MANTLE
    ggdof_elmt = reshape(ggdof_cm1(:,inode_elmt_cm1(:,i_elmt)),(/NEDOF1/))
    ggdof_elmt = ggdof_elmt-1 ! petsc index starts from 0
    ncount = 0; idof=-1; igdof=-1
    do i = 1,NEDOF1
      do j = 1,NEDOF1
      if (ggdof_elmt(i) >= 0.and.ggdof_elmt(j) >= 0.and.                          &
      storekmat_crust_mantle1(i,j,i_elmt) /= 0.0_kreal) then
        call MatSetValues(Amat1,1,ggdof_elmt(i),1,ggdof_elmt(j), &
        storekmat_crust_mantle1(i,j,i_elmt),ADD_VALUES,ierr)
        CHKERRA(ierr)
      endif
      enddo
    enddo
  enddo
  !do i_elmt=1,NSPEC_CRUST_MANTLE
  !  ggdof_elmt=reshape(ggdof_cm1(:,inode_elmt_cm1(:,i_elmt)),(/NEDOF1/))
  !  ggdof_elmt=ggdof_elmt-1 ! petsc index starts from 0
  !  ncount=0; idof=-1; igdof=-1
  !  do i=1,NEDOF1
  !    if (ggdof_elmt(i) >= 0) then
  !      ncount=ncount+1
  !      idof(ncount)=i
  !      igdof(ncount)=ggdof_elmt(i)
  !    endif
  !  enddo
  !  !if (myrank==0) write(*,*) 'hi homnath3in!',i_elmt,minval(igdof(1:ncount)), &
  !  !maxval(igdof(1:ncount)) !,storekmat_crust_mantle(idof(1:ncount),idof(1:ncount),i_elmt)
  !  call MatSetValues(Amat1,ncount,igdof(1:ncount),ncount,igdof(1:ncount), &
  !  storekmat_crust_mantle1(idof(1:ncount),idof(1:ncount),i_elmt),ADD_VALUES,ierr)
  !  CHKERRA(ierr)
  !enddo
  deallocate(storekmat_crust_mantle1)
  !if (myrank==0) write(*,*) 'CM kmat done1!'; call sync_all
  ! trinfinite
  do i_elmt = 1,NSPEC_TRINFINITE
    ggdof_elmt = reshape(ggdof_trinf1(:,inode_elmt_trinf1(:,i_elmt)),(/NEDOF1/))
    ggdof_elmt = ggdof_elmt-1 ! petsc index starts from 0
    ncount = 0; idof=-1; igdof=-1
    do i = 1,NEDOF1
      do j = 1,NEDOF1
      if (ggdof_elmt(i) >= 0.and.ggdof_elmt(j) >= 0.and.                          &
      storekmat_trinfinite1(i,j,i_elmt) /= 0.0_kreal) then
        call MatSetValues(Amat1,1,ggdof_elmt(i),1,ggdof_elmt(j), &
        storekmat_trinfinite1(i,j,i_elmt),ADD_VALUES,ierr)
        CHKERRA(ierr)
      endif
      enddo
    enddo
  enddo
  !do i_elmt=1,NSPEC_TRINFINITE
  !  ggdof_elmt=reshape(ggdof_trinf1(:,inode_elmt_trinf1(:,i_elmt)),(/NEDOF1/))
  !  ggdof_elmt=ggdof_elmt-1 ! petsc index starts from 0
  !  ncount=0; idof=-1; igdof=-1
  !  do i=1,NEDOF1
  !    if (ggdof_elmt(i) >= 0) then
  !      ncount=ncount+1
  !      idof(ncount)=i
  !      igdof(ncount)=ggdof_elmt(i)
  !    endif
  !  enddo
  !  call MatSetValues(Amat1,ncount,igdof(1:ncount),ncount,igdof(1:ncount), &
  !  storekmat_trinfinite1(idof(1:ncount),idof(1:ncount),i_elmt),ADD_VALUES,ierr)
  !  CHKERRA(ierr)
  !enddo
  deallocate(storekmat_trinfinite1)
  !if (myrank==0) write(*,*) 'TRINF kmat done1!'; call sync_all
  ! infinite
  do i_elmt = 1,NSPEC_INFINITE
    ggdof_elmt = reshape(ggdof_inf1(:,inode_elmt_inf1(:,i_elmt)),(/NEDOF1/))
    ggdof_elmt = ggdof_elmt-1 ! petsc index starts from 0
    ncount = 0; idof=-1; igdof=-1
    do i = 1,NEDOF1
      do j = 1,NEDOF1
      if (ggdof_elmt(i) >= 0.and.ggdof_elmt(j) >= 0.and.                          &
      storekmat_infinite1(i,j,i_elmt) /= 0.0_kreal) then
        call MatSetValues(Amat1,1,ggdof_elmt(i),1,ggdof_elmt(j), &
        storekmat_infinite1(i,j,i_elmt),ADD_VALUES,ierr)
        CHKERRA(ierr)
        !if (myrank==0) write(*,*) 'hello in INF1:',i_elmt,'/',NSPEC_INFINITE,i,j,NGLLCUBE_INF,NEDOF1!ggdof_elmt(i),ggdof_elmt(j)
      endif
      enddo
    enddo
  !  do i=1,NEDOF1
  !    if (ggdof_elmt(i) >= 0) then
  !      ncount=ncount+1
  !      idof(ncount)=i
  !      igdof(ncount)=ggdof_elmt(i)
  !    endif
  !  enddo
  !  call MatSetValues(Amat1,ncount,igdof(1:ncount),ncount,igdof(1:ncount), &
  !  storekmat_infinite1(idof(1:ncount),idof(1:ncount),i_elmt),ADD_VALUES,ierr)
  !  CHKERRA(ierr)
  enddo
  call sync_all
  !if (myrank==0) write(*,*) 'INF kmat done1:0!'; call sync_all
  deallocate(storekmat_infinite1)
  !if (myrank==0) write(*,*) 'INF kmat done1!'; call sync_all

  call MatAssemblyBegin(Amat1,MAT_FINAL_ASSEMBLY,ierr)
  CHKERRA(ierr)
  call MatAssemblyEnd(Amat1,MAT_FINAL_ASSEMBLY,ierr)
  CHKERRA(ierr)
  call MatSetOption(Amat1,MAT_SYMMETRIC,PETSC_TRUE,ierr);
  CHKERRA(ierr)
  !if (myrank==0) write(*,*) 'matrix setting & assembly complete11!'
  call sync_all
  call MatGetOwnershipRange(Amat1,istart,iend,ierr)
  CHKERRA(ierr)

  !! check diagonal of the matrix
  !call MatCreateVecs(Amat1,Adiag1,PETSC_NULL_OBJECT,ierr)
  !call MatGetDiagonal(Amat1,Adiag1,ierr)
  !
  !call VecCreateSeq(PETSC_COMM_SELF,neq1,lAdiag1,ierr)
  !
  !call VecScatterBegin(vscat1,Adiag1,lAdiag1,INSERT_VALUES,SCATTER_FORWARD,ierr)
  !CHKERRA(ierr)
  !call VecScatterEnd(vscat1,Adiag1,lAdiag1,INSERT_VALUES,SCATTER_FORWARD,ierr)
  !CHKERRA(ierr)
  !call VecGetSize(lAdiag1,n,ierr)
  !call VecGetArrayF90(lAdiag1,arrayAdiag1,ierr)
  !CHKERRA(ierr)
  !maxdiag1=maxscal(maxval(arrayAdiag1))
  !mindiag1=minscal(minval(arrayAdiag1))
  !write(*,*) 'maxdiag1 in petsc:',mindiag1,maxdiag1
  !call VecRestoreArrayF90(lAdiag1,arrayAdiag1,ierr)
  !CHKERRA(ierr)


  !allocate(cols(nzeros_max),vals(nzeros_max))
  allocate(cols(nzeros_max))
  !call MatGetRow(Amat1,istart,ncols,cols,vals,ierr);
  write(char_myrank,'(i4)')myrank
  outf_name='tmp/nonzeros'//trim(adjustl(char_myrank))
  open(1,file=outf_name,action='write',status='replace')
  do i = istart,iend-1
    cols=-1
    !call MatGetRow(Amat1,i,ncols,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr);
    call MatGetRow(Amat1,i,ncols,cols,PETSC_NULL_SCALAR,ierr);
    CHKERRA(ierr)
    !ncols_val=ncols
    !if (myrank==0) then
    !write(*,*) 'nzeros in 0th row:',myrank,istart,ncols
    ndiag=count(cols >= ig0.and.cols <= ig1)
    noffdiag = ncols-ndiag
    write(1,*)ndiag,noffdiag,ncols
    !endif
    !call MatRestoreRow(Amat1,istart,ncols,cols,vals,ierr);
    !call sync_all
    call MatRestoreRow(Amat1,i,ncols,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr);
    CHKERRA(ierr)
  enddo
  close(1)
  call sync_all()

  end subroutine petsc_set_matrix1

!
!===============================================================================
!

  subroutine petsc_set_vector1(rload1)

  use specfem_par, only: l2gdof1
  use specfem_par, only: IFLAG_IN_FICTITIOUS_CUBE,NSPEC_INNER_CORE, &
  NSPEC_OUTER_CORE,NSPEC_CRUST_MANTLE,NSPEC_TRINFINITE,NSPEC_INFINITE
  use specfem_par, only: NEDOF
  use specfem_par, only: NGLLCUBE_INF,NEDOF1
  use specfem_par_innercore, only: ggdof_ic1,storekmat_inner_core1, &
  idoubling_inner_core,inode_elmt_ic1
  use specfem_par_outercore, only: ggdof_oc1,storekmat_outer_core1, &
  inode_elmt_oc1
  use specfem_par_crustmantle, only: ggdof_cm1,storekmat_crust_mantle1, &
  inode_elmt_cm1
  use specfem_par_trinfinite, only: ggdof_trinf1,storekmat_trinfinite1, &
  inode_elmt_trinf1
  use specfem_par_infinite, only: ggdof_inf1,storekmat_infinite1,inode_elmt_inf1
  implicit none
  PetscScalar,intent(in) :: rload1(0:)
  PetscScalar      zero

  zero = 0.0
  call VecSet(bvec1,zero,ierr)
  call VecSetValues(bvec1,neq1,l2gdof1(1:),rload1(1:),ADD_VALUES,ierr);

  ! assemble vector
  call VecAssemblyBegin(bvec1,ierr)
  call VecAssemblyEnd(bvec1,ierr)
  !if (myrank==0) write(*,*) 'vector setting & assembly complete!'

  end subroutine petsc_set_vector1

!
!===============================================================================
!

  subroutine petsc_set_backward_vector1(b_rload1)

  use specfem_par, only: l2gdof1, neq1
  implicit none
  PetscScalar,intent(in) :: b_rload1(0:)
  PetscScalar      zero

  zero = 0.0
  call VecSet(b_bvec1,zero,ierr)
  call VecSetValues(b_bvec1,neq1,l2gdof1(1:),b_rload1(1:),ADD_VALUES,ierr);

  ! assemble vector
  call VecAssemblyBegin(b_bvec1,ierr)
  call VecAssemblyEnd(b_bvec1,ierr)

  end subroutine petsc_set_backward_vector1

!
!===============================================================================
!

  subroutine petsc_solve1(sdata1,iter,ireason)

  implicit none
  PetscInt iter
  PetscScalar sdata1(:)

  PetscInt ireason


  call KSPSolve(ksp1,bvec1,xvec1,ierr)
  call KSPGetConvergedReason(ksp1,ireason,ierr)
  call KSPGetIterationNumber(ksp1,iter,ierr)

  call scatter_globalvec1(xvec1, sdata1)

  end subroutine petsc_solve1

!
!===============================================================================
!

  subroutine petsc_backward_solve1(b_sdata1,b_iter,b_ireason)

  implicit none
  PetscInt b_iter
  PetscScalar b_sdata1(:)
  PetscInt b_ireason

  call KSPSolve(b_ksp1,b_bvec1,b_xvec1,ierr)
  call KSPGetConvergedReason(b_ksp1,b_ireason,ierr)
  call KSPGetIterationNumber(b_ksp1,b_iter,ierr)

  call scatter_globalvec1_backward(b_xvec1, b_sdata1)

  end subroutine petsc_backward_solve1

!
!===============================================================================
!

!  subroutine scatter_globalvec1(global_vec, larray, l_vec1, l_vscat1)
!  ! l_vec1 is local_vec1 or b_local_vec1 for forward/adjoint solver
!  ! l_vscat1 is vscat1 or b_vscat1
!  implicit none
!  ! I/O variables
!  VecScatter l_vscat1
!  Vec l_vec1
!  Vec global_vec
!  PetscScalar larray(:)
!  ! Local variables
!  PetscInt n
!  PetscScalar,pointer :: array_data(:)
!
!  call VecScatterBegin(l_vscat1,global_vec,l_vec1,INSERT_VALUES,SCATTER_FORWARD,ierr)
!  CHKERRA(ierr)
!  call VecScatterEnd(l_vscat1,global_vec,l_vec1,INSERT_VALUES,SCATTER_FORWARD,ierr)
!  CHKERRA(ierr)
!  call VecGetSize(l_vec1,n,ierr)
!  call VecGetArrayF90(l_vec1,array_data,ierr)
!  CHKERRA(ierr)
!  larray(1:n)=array_data(1:n)
!  call VecRestoreArrayF90(l_vec1,array_data,ierr)
!  CHKERRA(ierr)
!  return
!  end subroutine scatter_globalvec1

!
!===============================================================================
!

  subroutine scatter_globalvec1(global_vec,larray)

  implicit none

  Vec global_vec
  PetscScalar larray(:)
  PetscInt n

  PetscScalar,pointer :: array_data(:)

  call VecScatterBegin(vscat1,global_vec,local_vec1,INSERT_VALUES,SCATTER_FORWARD,ierr)
  CHKERRA(ierr)
  call VecScatterEnd(vscat1,global_vec,local_vec1,INSERT_VALUES,SCATTER_FORWARD,ierr)
  CHKERRA(ierr)
  call VecGetSize(local_vec1,n,ierr)
  call VecGetArrayF90(local_vec1,array_data,ierr)
  CHKERRA(ierr)
  larray(1:n)=array_data(1:n)
  call VecRestoreArrayF90(local_vec1,array_data,ierr)
  CHKERRA(ierr)
  return

  end subroutine scatter_globalvec1

!
!===============================================================================
!

  subroutine scatter_globalvec1_backward(b_global_vec,b_larray)

  implicit none

  Vec b_global_vec
  PetscScalar b_larray(:)
  PetscInt b_n

  PetscScalar,pointer :: b_array_data(:)

  call VecScatterBegin(b_vscat1, b_global_vec, b_local_vec1,INSERT_VALUES,SCATTER_FORWARD,ierr)
  CHKERRA(ierr)
  call VecScatterEnd(b_vscat1, b_global_vec, b_local_vec1,INSERT_VALUES,SCATTER_FORWARD,ierr)
  CHKERRA(ierr)
  call VecGetSize(b_local_vec1, b_n,ierr)
  call VecGetArrayF90(b_local_vec1, b_array_data,ierr)
  CHKERRA(ierr)
  b_larray(1:b_n)=b_array_data(1:b_n)
  call VecRestoreArrayF90(b_local_vec1, b_array_data,ierr)
  CHKERRA(ierr)
  return

  end subroutine scatter_globalvec1_backward

!
!===============================================================================
!

  subroutine petsc_zero_initialguess1()

  implicit none
  PetscScalar      zero

  zero = 0.0
  call VecSet(xvec1,zero,ierr)

  ! assemble vector
  call VecAssemblyBegin(xvec1,ierr)
  call VecAssemblyEnd(xvec1,ierr)

  end subroutine petsc_zero_initialguess1

!
!===============================================================================
!

  subroutine petsc_zero_backwards_initialguess1()

  implicit none
  PetscScalar      zero

  zero = 0.0
  call VecSet(b_xvec1,zero,ierr)

  ! assemble vector
  call VecAssemblyBegin(b_xvec1,ierr)
  call VecAssemblyEnd(b_xvec1,ierr)

  end subroutine petsc_zero_backwards_initialguess1

!
!===============================================================================
!

!  subroutine petsc_set_initialguess1(loc_pgrav1)
!  implicit none
!  PetscScalar,intent(in) :: loc_pgrav1(0:)
!  PetscScalar      zero
!
!  zero=0.0
!  call VecSet(bvec1,zero,ierr)
!  call VecSetValues(bvec1,neq1,l2gdof1(1:),loc_pgrav1(1:),INSERT_VALUES,ierr);
!
!  ! assemble vector
!  call VecAssemblyBegin(bvec1,ierr)
!  call VecAssemblyEnd(bvec1,ierr)
!
!  end subroutine petsc_set_initialguess1

!
!===============================================================================
!

!  subroutine petsc_set_backward_initialguess1(loc_pgrav1)
!    implicit none
!    PetscScalar,intent(in) :: loc_pgrav1(0:)
!    PetscScalar      zero
!
!    zero=0.0
!    call VecSet(b_bvec1,zero,ierr)
!    call VecSetValues(b_bvec1, neq1, l2gdof1(1:), loc_pgrav1(1:), INSERT_VALUES, ierr);
!
!    ! assemble vector
!    call VecAssemblyBegin(b_bvec1,ierr)
!    call VecAssemblyEnd(b_bvec1,ierr)
!
!  end subroutine petsc_set_backward_initialguess1

!
!===============================================================================
!

  subroutine petsc_finalize1()

  use specfem_par, only: SIMULATION_TYPE

  implicit none

  ! Free work space.  All PETSc objects should be destroyed when they
  ! are no longer needed.

  call VecDestroy(xvec1,ierr)
  call VecDestroy(uvec1,ierr)
  call VecDestroy(bvec1,ierr)
  call MatDestroy(Amat1,ierr)
  call KSPDestroy(ksp1,ierr)
  call VecScatterDestroy(vscat1,ierr)
  call PetscFinalize(ierr)

  if (SIMULATION_TYPE == 3) then
    call VecDestroy(b_xvec1,ierr)
    call VecDestroy(b_bvec1,ierr)
    call KSPDestroy(b_ksp1,ierr)
  call VecScatterDestroy(b_vscat1,ierr)
  endif

  end subroutine petsc_finalize1


!===============================================================================
! Level-2 solver
!===============================================================================

  subroutine petsc_initialize()

  implicit none
  PetscInt :: istart,iend
  PetscInt :: nzeros_max,nzeros_min
  PetscInt, allocatable :: nzeros(:)
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Beginning of program
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  !call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-n',ngdof,flg,ierr)
  !if (myrank==0) write(*,*) 'hi0!'
  !call sync_all
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Compute the matrix and right-hand-side vector that define
  ! the linear system, Ax = b.
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ! Create matrix. When using MatCreate(), the matrix format can
  ! be specified at runtime.

  ! count number of nonzeros per row
  allocate(nzeros(neq))
  nzeros = 0
  nzeros(krow_sparse)=nzeros(krow_sparse)+1
  nzeros_max=maxvec(nzeros)
  nzeros_min=minvec(nzeros)
  !nzeros_max=2*nzeros_max
  !nzeros=nzeros
  if (myrank == 0) write(*,*) 'ngdof:',ngdof,' nzeros_max:',nzeros_max,' nzeros_min:',nzeros_min
  call MatCreate(PETSC_COMM_WORLD,Amat,ierr)
  call MatSetType(Amat,MATMPIAIJ,ierr)
  CHKERRA(ierr)
  call MatSetSizes(Amat,PETSC_DECIDE,PETSC_DECIDE,ngdof,ngdof,ierr)
  CHKERRA(ierr)
  ! preallocation
  !call MatMPIAIJSetPreallocation(Amat,nzeros_max,PETSC_NULL_INTEGER,nzeros_max, &
  !PETSC_NULL_INTEGER,ierr)
  call MatMPIAIJSetPreallocation(Amat,nzeros_max,nzeros,nzeros_max, &
  20*nzeros,ierr)
  CHKERRA(ierr)
  !call MatSeqAIJSetPreallocation(Amat,nzeros_max,nzeros,ierr)
  !CHKERRA(ierr)

  !call MatCreateMPIAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,ngdof,ngdof, &
  !nzeros_max,PETSC_NULL_INTEGER,nzeros_max,PETSC_NULL_INTEGER,Amat,ierr)

  !if (myrank==0) write(*,*) 'Matrix size:',size(Amat,1),size(Amat,2)
  call MatSetFromOptions(Amat,ierr)
  CHKERRA(ierr)
  !call MatSetUp(Amat,ierr)
  !if (myrank==0) write(*,*) 'ierr1:',ierr

  call MatGetOwnershipRange(Amat,istart,iend,ierr)
  CHKERRA(ierr)

  !if (myrank==0) write(*,*) 'ierr2:',ierr,istart,iend,iend-istart,minval(nzeros),maxval(nzeros)
  if (myrank == 0) write(*,*) 'actual global index range:',minval(kgrow_sparse),maxval(kgrow_sparse)
  write(*,*) 'global index:',myrank,istart,iend,iend-istart
  call sync_all
  !if (myrank==0) write(*,*) 'ierr3:',(iend-istart)*nzeros
  !if (myrank==0) write(*,*) 'ierr4:',iend,sum(nzeros),sum((iend-istart)*nzeros)
  deallocate(nzeros)
  !call sync_all
  call sync_all
  if (myrank == 0) write(*,*) 'matrix'
  !call sync_all


  ! Create vectors.  Note that we form 1 vector from scratch and
  ! then duplicate as needed.

  !call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,ngdof,xvec,ierr)
  call VecCreate(PETSC_COMM_WORLD,xvec,ierr)
  call VecSetSizes(xvec,PETSC_DECIDE,ngdof,ierr)
  call VecSetFromOptions(xvec,ierr)
  call VecDuplicate(xvec,bvec,ierr)
  call VecDuplicate(xvec,uvec,ierr)
  if (myrank == 0) write(*,*) 'vector'

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Create the linear solver and set various options
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ! Create linear solver context

  call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
  ! Set operators. Here the matrix that defines the linear system
  ! also serves as the preconditioning matrix.
  !call KSPSetOperators(ksp,Amat,Amat,SAME_PRECONDITIONER,ierr) ! version < 3.5
  call KSPSetOperators(ksp,Amat,Amat,ierr) ! version >= 3.5

  call KSPSetType(ksp,KSPCG,ierr);
  if (myrank == 0) write(*,*) 'ksp0'
  call KSPGetPC(ksp,pc,ierr)
  call PCSetType(pc,PCHYPRE,ierr)
  if (myrank == 0) write(*,*) 'ksp1'
  call KSPSetTolerances(ksp,KSP_RTOL,KSP_ATOL,KSP_DTOL,KSP_MAXITER,ierr)
  CHKERRA(ierr)
  if (myrank == 0) write(*,*) 'ksp2'

  !  Set runtime options, e.g.,
  !    -ksp_type < type> -pc_type < type> -ksp_monitor -ksp_KSP_RTOL < KSP_RTOL>
  !  These options will override those specified above as long as
  !  KSPSetFromOptions() is called _after_ any other customization
  !  routines.
  call KSPSetFromOptions(ksp,ierr)

  end subroutine petsc_initialize

!
!===============================================================================
!

  subroutine petsc_set_matrix()

  use specfem_par, only: NEDOF,IFLAG_IN_FICTITIOUS_CUBE,NSPEC_INNER_CORE, &
  NSPEC_OUTER_CORE,NSPEC_CRUST_MANTLE,NSPEC_TRINFINITE,NSPEC_INFINITE
  use specfem_par_innercore, only: ggdof_ic,storekmat_inner_core, &
  idoubling_inner_core,inode_elmt_ic
  use specfem_par_outercore, only: ggdof_oc,storekmat_outer_core,inode_elmt_oc
  use specfem_par_crustmantle, only: ggdof_cm,storekmat_crust_mantle,inode_elmt_cm
  use specfem_par_trinfinite, only: ggdof_trinf,storekmat_trinfinite,inode_elmt_trinf
  use specfem_par_infinite, only: ggdof_inf,storekmat_infinite,inode_elmt_inf

  implicit none
  integer :: i,i_elmt,j,ncount
  integer :: ggdof_elmt(NEDOF),idof(NEDOF),igdof(NEDOF)
  ! Set and assemble matrix.
  !  - Note that MatSetValues() uses 0-based row and column numbers
  !  in Fortran as well as in C (as set here in the array "col").
    ! stage 0: store all elements

  call MatZeroEntries(Amat,ierr)
  CHKERRA(ierr)

  ! inner core
  do i_elmt = 1,NSPEC_INNER_CORE
    if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE)cycle
    ggdof_elmt = reshape(ggdof_ic(:,inode_elmt_ic(:,i_elmt)),(/NEDOF/))
    ggdof_elmt = ggdof_elmt-1 ! petsc index starts from 0
    ncount = 0; idof=-1; igdof=-1
    do i = 1,NEDOF
      do j = 1,NEDOF
      if (ggdof_elmt(i) >= 0.and.ggdof_elmt(j) >= 0.and.                          &
      storekmat_inner_core(i,j,i_elmt) /= 0.0_kreal) then
        !ncount=ncount+1
        !idof(ncount)=i
        !igdof(ncount)=ggdof_elmt(i)
        call MatSetValues(Amat,1,ggdof_elmt(i),1,ggdof_elmt(j), &
        storekmat_inner_core(i,j,i_elmt),ADD_VALUES,ierr)
        CHKERRA(ierr)
      endif
      enddo
    enddo
    !call MatSetValues(Amat,ncount,igdof(1:ncount),ncount,igdof(1:ncount), &
    !storekmat_inner_core(idof(1:ncount),idof(1:ncount),i_elmt),ADD_VALUES,ierr)
    !CHKERRA(ierr)
  enddo
  !! inner core
  !do i_elmt=1,NSPEC_INNER_CORE
  !  if (idoubling_inner_core(i_elmt)==IFLAG_IN_FICTITIOUS_CUBE)cycle
  !  ggdof_elmt=reshape(ggdof_ic(:,inode_elmt_ic(:,i_elmt)),(/NEDOF/))
  !  ggdof_elmt=ggdof_elmt-1 ! petsc index starts from 0
  !  ncount=0; idof=-1; igdof=-1
  !  do i=1,NEDOF
  !    if (ggdof_elmt(i) >= 0) then
  !      ncount=ncount+1
  !      idof(ncount)=i
  !      igdof(ncount)=ggdof_elmt(i)
  !    endif
  !  enddo
  !  call MatSetValues(Amat,ncount,igdof(1:ncount),ncount,igdof(1:ncount), &
  !  storekmat_inner_core(idof(1:ncount),idof(1:ncount),i_elmt),ADD_VALUES,ierr)
  !  CHKERRA(ierr)
  !enddo
  deallocate(storekmat_inner_core)
  !if (myrank==0) write(*,*) 'IC kmat done!'; call sync_all
  ! outer core
  do i_elmt = 1,NSPEC_OUTER_CORE
    ggdof_elmt = reshape(ggdof_oc(:,inode_elmt_oc(:,i_elmt)),(/NEDOF/))
    ggdof_elmt = ggdof_elmt-1 ! petsc index starts from 0
    ncount = 0; idof=-1; igdof=-1
    do i = 1,NEDOF
      do j = 1,NEDOF
      if (ggdof_elmt(i) >= 0.and.ggdof_elmt(j) >= 0.and.                          &
      storekmat_outer_core(i,j,i_elmt) /= 0.0_kreal) then
        if (myrank == 0) write(*,*) 'hello in OC:',i_elmt,ggdof_elmt(i),ggdof_elmt(j)
        call sync_all
        call MatSetValues(Amat,1,ggdof_elmt(i),1,ggdof_elmt(j), &
        storekmat_outer_core(i,j,i_elmt),ADD_VALUES,ierr)
        CHKERRA(ierr)
      endif
      enddo
    enddo
  enddo
  !! outer core
  !do i_elmt=1,NSPEC_OUTER_CORE
  !  ggdof_elmt=reshape(ggdof_oc(:,inode_elmt_oc(:,i_elmt)),(/NEDOF/))
  !  ggdof_elmt=ggdof_elmt-1 ! petsc index starts from 0
  !  ncount=0; idof=-1; igdof=-1
  !  do i=1,NEDOF
  !    if (ggdof_elmt(i) >= 0) then
  !      ncount=ncount+1
  !      idof(ncount)=i
  !      igdof(ncount)=ggdof_elmt(i)
  !    endif
  !  enddo
  !  call MatSetValues(Amat,ncount,igdof(1:ncount),ncount,igdof(1:ncount), &
  !  storekmat_outer_core(idof(1:ncount),idof(1:ncount),i_elmt),ADD_VALUES,ierr)
  !  CHKERRA(ierr)
  !enddo
  deallocate(storekmat_outer_core)
  !if (myrank==0) write(*,*) 'OC kmat done!'; call sync_all
  ! crust mantle
  do i_elmt = 1,NSPEC_CRUST_MANTLE
    ggdof_elmt = reshape(ggdof_cm(:,inode_elmt_cm(:,i_elmt)),(/NEDOF/))
    ggdof_elmt = ggdof_elmt-1 ! petsc index starts from 0
    ncount = 0; idof=-1; igdof=-1
    do i = 1,NEDOF
      if (ggdof_elmt(i) >= 0) then
        ncount = ncount+1
        idof(ncount)=i
        igdof(ncount)=ggdof_elmt(i)
      endif
    enddo
    !if (myrank==0) write(*,*) 'hi homnath3in!',i_elmt,minval(igdof(1:ncount)), &
    !maxval(igdof(1:ncount)) !,storekmat_crust_mantle(idof(1:ncount),idof(1:ncount),i_elmt)
    call MatSetValues(Amat,ncount,igdof(1:ncount),ncount,igdof(1:ncount), &
    storekmat_crust_mantle(idof(1:ncount),idof(1:ncount),i_elmt),ADD_VALUES,ierr)
    CHKERRA(ierr)
  enddo
  deallocate(storekmat_crust_mantle)
  !if (myrank==0) write(*,*) 'CM kmat done!'; call sync_all
  ! trinfinite
  do i_elmt = 1,NSPEC_TRINFINITE
    ggdof_elmt = reshape(ggdof_trinf(:,inode_elmt_trinf(:,i_elmt)),(/NEDOF/))
    ggdof_elmt = ggdof_elmt-1 ! petsc index starts from 0
    ncount = 0; idof=-1; igdof=-1
    do i = 1,NEDOF
      if (ggdof_elmt(i) >= 0) then
        ncount = ncount+1
        idof(ncount)=i
        igdof(ncount)=ggdof_elmt(i)
      endif
    enddo
    call MatSetValues(Amat,ncount,igdof(1:ncount),ncount,igdof(1:ncount), &
    storekmat_trinfinite(idof(1:ncount),idof(1:ncount),i_elmt),ADD_VALUES,ierr)
    CHKERRA(ierr)
  enddo
  deallocate(storekmat_trinfinite)
  !if (myrank==0) write(*,*) 'TRINF kmat done!'; call sync_all
  ! infinite
  do i_elmt = 1,NSPEC_INFINITE
    ggdof_elmt = reshape(ggdof_inf(:,inode_elmt_inf(:,i_elmt)),(/NEDOF/))
    ggdof_elmt = ggdof_elmt-1 ! petsc index starts from 0
    ncount = 0; idof=-1; igdof=-1
    do i = 1,NEDOF
      if (ggdof_elmt(i) >= 0) then
        ncount = ncount+1
        idof(ncount)=i
        igdof(ncount)=ggdof_elmt(i)
      endif
    enddo
    call MatSetValues(Amat,ncount,igdof(1:ncount),ncount,igdof(1:ncount), &
    storekmat_infinite(idof(1:ncount),idof(1:ncount),i_elmt),ADD_VALUES,ierr)
    CHKERRA(ierr)
  enddo
  deallocate(storekmat_infinite)
  !if (myrank==0) write(*,*) 'INF kmat done!'; call sync_all

  call MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY,ierr)
  !if (myrank==0) write(*,*) 'matrix setting & assembly complete!'; call sync_all

  end subroutine petsc_set_matrix

!
!===============================================================================
!

  subroutine petsc_set_vector()

  use specfem_par, only: l2gdof,load
  implicit none
  PetscScalar      zero !,none,one
  ! Set exact solution; then compute right-hand-side vector.
  !none=-1.0
  !one=1.0
  zero = 0.0
  call VecSet(bvec,zero,ierr)
  call VecSetValues(bvec,neq,l2gdof(1:),load(1:),ADD_VALUES,ierr);

  ! assemble vector
  call VecAssemblyBegin(bvec,ierr)
  call VecAssemblyEnd(bvec,ierr)
  !if (myrank==0) write(*,*) 'vector setting & assembly complete!'

  end subroutine petsc_set_vector

!
!===============================================================================
!

  subroutine petsc_solve(sdata,cg_iter)

  implicit none
  PetscScalar sdata(:)
  PetscInt    cg_iter
  PetscInt    ireason

  call KSPSolve(ksp,bvec,xvec,ierr)

  ! View solver info; we could instead use the option -ksp_view
  call KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr)

  !-------------------------------------------------------------------------------
  ! Check solution and clean up
  !-------------------------------------------------------------------------------
  ! Check the error
  !call VecAXPY(xvec,none,uvec,ierr)
  !call VecNorm(xvec,NORM_2,norm,ierr)
  call KSPGetConvergedReason(ksp,ireason,ierr)
  call KSPGetIterationNumber(ksp,cg_iter,ierr)
  if (myrank < 1) then
    write(*,*) 'converged reason',ireason
    write(*,*) 'Iterations:',cg_iter
  endif
  !if (norm > 1.e-12) then
  !  write(*,'(a,e11.4,a,i5)')'Norm of error:',norm,', Iterations:',its
  !else
  !  write(*,'(a,i5,a)')'Norm of error < 1.e-12, Iterations:',its
  !endif

  end subroutine petsc_solve

!
!===============================================================================
!

  subroutine petsc_finalize()

  implicit none

  ! Free work space.  All PETSc objects should be destroyed when they
  ! are no longer needed.

  call VecDestroy(xvec,ierr)
  call VecDestroy(uvec,ierr)
  call VecDestroy(bvec,ierr)
  call MatDestroy(Amat,ierr)
  call KSPDestroy(ksp,ierr)
  call PetscFinalize(ierr)

  end subroutine petsc_finalize

end module solver_petsc

#endif

