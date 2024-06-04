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

!AUTHORS:
!Hom Nath Gharti
!Stefano Zhampini
!REFERENCE:
!PETSC documentation

! 1. note on PETSc macros:
!
! the PETSc base type file petscsys.h defines an error checking macro CHKERRA(ierr) as
! long free-from version:
!#define CHKERRA(ierr) if (ierr /= 0) then;call PetscErrorF(ierr,__LINE__,__FILE__);\
!                      call MPIU_Abort(PETSC_COMM_SELF,ierr);endif
! short version:
!#define CHKERRA(ierr) if (ierr /= 0) then;call PetscErrorF(ierr);\
!                      call MPIU_Abort(PETSC_COMM_SELF,ierr);endif
!
! with gfortran, it uses the long free-form version that leads to compilation errors like
!   "Error: Line truncated at (1) [-Werror=line-truncation]"
! whenever this macro is indented, for example in an if-statement like:
! if (..) then
!    ..
!    CHKERRA(ierr)
! endif
! it only works if the macro gets indented by 2 spaces, but not more :(
!
! PETSc suggests to call its functions with the PetscCall macro, for example
!   PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
! where the corresponding macro is defined as:
!   #define PetscCallA(func) call func; CHKERRA(ierr)
! unfortunately, this leads to the same compilation errors with exceeding line length limits.
!
! with the Fortran preprocessors, macro definitions unfortunately can't span over multiple lines
! and thus are limited by length.
! thus, we will define an error checking subroutine instead in this module,
! and call it by our own marco to work-around this limit.

#define CHECK_PETSC_ERROR(ierr) call check_petsc_err(ierr,__LINE__)


! 2. note on PETSc types and function calls:
!
! the Petsc types might differ a bit on different systems and compilers.
! for example, calling the function
!   VecSetValues(Vec x, PetscInt ni, const PetscInt ix[], const PetscScalar y[], InsertMode iora)
! with a real(kind=CUSTOM_REAL) gravload(:) array for the PetscScalar array y[] will cause a compilation error on macOS like:
!   call VecSetValues(bvec,neq,l2gdof(1:),gravload(1:),ADD_VALUES,ierr)
!   ..
!   "Error: There is no specific subroutine for the generic 'vecsetvalues' at (1)"
!
! we thus use copy-arrays and explicit type conversions to work around this:
!   PetsScalar :: y(size(gravload(1:))
!   y(:) = gravload(1:)
!   call VecSetValues(bvec,neq,l2gdof(1:),y,ADD_VALUES,ierr)
! well, until a better solution is found.


module siem_solver_petsc

#ifdef USE_PETSC
! PETSc
! all of PETSc
#include "petsc/finclude/petsc.h"
  use petsc
! base types
#include "petsc/finclude/petscsys.h"
  use petscsys
! Vec package
#include "petsc/finclude/petscvec.h"
  use petscvec
! Mat package
#include "petsc/finclude/petscmat.h"
  use petscmat
! IS (index set) package
#include "petsc/finclude/petscis.h"
  use petscis
! Krylov subspace method package
#include "petsc/finclude/petscksp.h"
  use petscksp

! preconditioner package
!#include "petsc/finclude/petscpc.h"


! PETSc modules
  !use petsc, only: PETSC_COMM_WORLD,PETSC_COMM_SELF,PETSC_VIEWER_STDOUT_WORLD,PETSC_COPY_VALUES, &
  !                 PETSC_NULL_CHARACTER,PETSC_NULL_SCALAR,PETSC_NULL_OPTIONS,PETSC_NULL_INTEGER, &
  !                 PETSC_DECIDE,PETSC_TRUE,PETSC_FALSE
  !use petscmat, only: tMat,MAT_FINAL_ASSEMBLY,MAT_SYMMETRIC,MAT_SHIFT_POSITIVE_DEFINITE,MAT_SPD
  !use petscvec, only: tVec,tVecScatter,ADD_VALUES,INSERT_VALUES,SCATTER_FORWARD
  !use petscksp, only: tKSP
  !use petscpc, only: tPC
  !use petscis, only: tIS

#endif

  use constants, only: myrank,IMAIN,CUSTOM_REAL

  use constants_solver, only: &
    nproc => NPROCTOT_VAL

  use specfem_par_full_gravity, only: &
    neq,ngdof,nsparse,krow_sparse,kcol_sparse,kgrow_sparse,kgcol_sparse, &
    neq1,ngdof1,nsparse1,krow_sparse1,kcol_sparse1,kgrow_sparse1,kgcol_sparse1, &
    l2gdof1

  use specfem_par_full_gravity, only: &
    num_interfaces_inner_core1, &
    max_nibool_interfaces_inner_core1,my_neighbors_inner_core1, &
    nibool_interfaces_inner_core1,ibool_interfaces_inner_core1, &
    num_interfaces_outer_core1, &
    max_nibool_interfaces_outer_core1,my_neighbors_outer_core1, &
    nibool_interfaces_outer_core1,ibool_interfaces_outer_core1, &
    num_interfaces_crust_mantle1, &
    max_nibool_interfaces_crust_mantle1,my_neighbors_crust_mantle1, &
    nibool_interfaces_crust_mantle1,ibool_interfaces_crust_mantle1, &
    num_interfaces_trinfinite1, &
    max_nibool_interfaces_trinfinite1,my_neighbors_trinfinite1, &
    nibool_interfaces_trinfinite1,ibool_interfaces_trinfinite1, &
    num_interfaces_infinite1, &
    max_nibool_interfaces_infinite1,my_neighbors_infinite1, &
    nibool_interfaces_infinite1,ibool_interfaces_infinite1

  use siem_math_library_mpi, only: maxvec,minvec

  implicit none

  private

#ifdef USE_PETSC
  ! to match PETSc types:
  !PetscInt :: i
  !PetscReal :: x
  !PetscBool :: b
  !PetscErrorCode :: ierr
  !
  !integer, parameter :: petsc_int = kind(i)        ! integer kind to be used in SIEM solver petsc
  !integer, parameter :: petsc_real = kind(x)       ! real kind to be used in SIEM solver petsc
  !integer, parameter :: petsc_bool = kind(b)       ! boolean kind to be used in SIEM solver petsc
  !integer, parameter :: petsc_err = kind(ierr)     ! error kind to be used in SIEM solver petsc
  !
  !integer(petsc_int), parameter :: COMMAND = 0
  !integer(petsc_bool) :: flg
  !integer(petsc_int) :: ival
  !real(petsc_real) :: val
  !
  ! or directly use PETSC type macros...

  ! Parameters
  ! Krylov subspace method (KSP) tolerances
  !   types required by KSPSetTolerances:
  !   KSPSetTolerances(KSP ksp, PetscReal rtol, PetscReal abstol, PetscReal dtol, PetscInt maxits)
  !
  ! Level-1 KSP solver
  PetscInt, parameter :: KSP_MAXITER1 = 3000
  PetscReal, parameter :: KSP_RTOL1 = 1.0e-7
  PetscReal, parameter :: KSP_ATOL1 = 1.0e-30
  PetscReal, parameter :: KSP_DTOL1 = 1.0e30

  ! Level-2 KSP solver
  PetscInt, parameter :: KSP_MAXITER = 3000
  PetscReal, parameter :: KSP_RTOL = 1.0e-7
  PetscReal, parameter :: KSP_ATOL = 1.0e-30
  PetscReal, parameter :: KSP_DTOL = 1.0e30

  ! solver type
  PetscInt, parameter     :: COMMAND = 0, CG = 1, SUPERLU = 2, MUMPS = 3

  ! module variables
  !PetscBool               :: flg,flg_ch,flg_lu,flg_ilu
  !PetscInt                :: ival,icntl
  !PetscReal               :: val

  ! Level-1 solver
  type(tVec)              :: xvec1,bvec1,uvec1,local_vec1
  type(tMat)              :: Amat1 !,Fmat1
  type(tKSP)              :: ksp1
  type(tPC)               :: pc1
  !PetscInt                :: solver_type1 ! solver type

  ! For communications from local to global
  type(tVecScatter)       :: vscat1 !,pscat1

  ! Stores l2g map info
  !ISLocalToGlobalMapping  :: l2gmap

  !PetscBool              ::flg

  ! ADJOINT SIMULATIONS
  ! FOR NOW WE ASSUME THAT THE FORWARD AND ADJOINT SIMULATIONS ARE SOLVED
  ! WITH A CG METHOD, AND THAT ONLY LEVEL 1 SOLVER IS USED
  !
  ! Adjoint Level-1 solver
  ! notes: no b_uvec1 since it seems to just be created and destroyed
  !        no bAmat1 since we can use the same stiffness matrix
  !        I think we can probably use b_pc1 = pc1 but unsure
  type(tVec)              :: b_xvec1, b_bvec1, b_local_vec1
  type(tKSP)              :: b_ksp1
  type(tPC)               :: b_pc1

  ! For communications from local to global
  type(tVecScatter)       :: b_vscat1

  ! Level-2 solver
  type(tVec)              :: xvec,bvec,uvec
  type(tMat)              :: Amat
  type(tKSP)              :: ksp
  type(tPC)               :: pc
  PetscErrorCode          :: ierr
  PetscInt                :: nzeros_max,nzeros_min
  PetscInt                :: ig0,ig1
#endif

  ! public function
  public :: petsc_initialize1
  public :: petsc_initialize
  public :: petsc_finalize1
  public :: petsc_finalize
  public :: petsc_set_vector1
  public :: petsc_set_vector
  public :: petsc_set_matrix1
  public :: petsc_set_matrix
  public :: petsc_solve1
  public :: petsc_solve
  public :: petsc_zero_initialguess1
  public :: petsc_zero_backwards_initialguess1

contains

!===============================================================================
! Level-1 solver
!===============================================================================
  subroutine petsc_initialize1()

#ifdef USE_PETSC
  use constants, only: NNDOF

  use specfem_par, only: ADD_TRINF,SIMULATION_TYPE

  use specfem_par_full_gravity, only: ggdof_ic1,ggdof_oc1,ggdof_cm1,ggdof_trinf1,ggdof_inf1

  implicit none
  type(tVec) :: nzeror_gvec1,nzeror_dvec1,nzeror_ovec1,iproc_gvec1, &
                interface_gvec1,ninterface_dvec1,ninterface_ovec1,nself_gvec1
  PetscInt :: i,istart,iend,n
  !PetscInt :: nzerosoff_max
  PetscInt :: inzeros_max,inzeros_min
  PetscInt,allocatable :: nzeros(:),ig_array1(:)
  PetscScalar,allocatable :: rproc_array1(:)
  PetscScalar :: rval
  type(tIS) :: global_is, local_is
  type(tIS) :: b_global_is, b_local_is

  PetscInt :: igdof,ind,maxrank0,ng,ng0,ng1,np0
  PetscInt,allocatable :: inzeror_array1(:),iproc_array1(:),nzeros_row(:)
  PetscInt,allocatable :: nnzero_diag1(:),nnzero_offdiag1(:)
  PetscScalar,pointer :: nzeror_array1(:)
  PetscScalar,pointer :: nzeror_darray1(:),nzeror_oarray1(:),rnself_array1(:)
  PetscReal :: fac_ni,max_ni,pmax,pmin,rnid,rnioffd,rnd,rnoffd

  PetscInt :: ir,ic,igr,igc,ir0,ic0,igr0,igc0
  PetscInt :: nd,noffd
  PetscInt :: i_bool,i_ndof

  PetscInt :: nibool
  PetscInt,allocatable :: ibool_interface(:),isg_interface(:),nself_array1(:)
  PetscInt,allocatable :: ninterface_darray1(:),ninterface_oarray1(:)
  PetscScalar,allocatable :: rg_interface(:),rnself_lgarray1(:)
  PetscScalar,pointer :: rninterface_darray1(:),rninterface_oarray1(:)

  ! memory info
  PetscLogDouble :: bytes

  !character(len=10) :: char_myrank
  character(len=128) :: version

  ! timing
  double precision :: tstart,tCPU
  double precision, external :: wtime

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    initializing PETSc Level-1 solver'
    call flush_IMAIN()
  endif

  ! timing
  tstart = wtime()

  if (myrank == 0) print *
  if (myrank == 0) print *,'PETSc solver: ---------- Initialise PETSC: ---------- '

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr); CHECK_PETSC_ERROR(ierr)

  ! version info
  call PetscGetVersion(version,ierr); CHECK_PETSC_ERROR(ierr)
  if (myrank == 0) print *,'PETSc solver: version: ',trim(version)

  ! count number of nonzeros per row
  allocate(nzeros(neq1))
  nzeros(:) = 0
  do i = 1,nsparse1
    nzeros(krow_sparse1(i)) = nzeros(krow_sparse1(i)) + 1
  enddo

  nzeros_max = maxvec(nzeros)
  nzeros_min = minvec(nzeros)

  !nzerosoff_max = nzeros_max
  !nzeros_max = 4*nzeros_max
  !nzeros = nzeros
  !nzeros = 5*nzeros

  if (myrank == 0) print *,'PETSc solver: nzeros in 1st index:',nzeros(1)
  if (myrank == 0) print *,'PETSc solver: ngdof1:',ngdof1,' nzeros_max:',nzeros_max,' nzeros_min:',nzeros_min, &
                                         ' count:',count(krow_sparse1 == 1)
  ! free temporary array
  deallocate(nzeros)

  ! precompute ownership range OR partion layout
  ng1 = ngdof1/nproc
  ng0 = ceiling(real(ngdof1)/real(nproc))

  np0 = ngdof1 - nproc*ng1

  if (np0 == 0) then
    ! ng0=ng1
    ! all processors have equal gdofs
    ng = ng0
    ig0 = myrank*ng0 ! 0-based index
    ig1 = ig0+ng0-1
  else if (np0 > 0) then
    ! first np0 processors have ng0 gdofs each and remaining processors have ng1
    ! gdofs each
    maxrank0 = np0-1 ! myrank is 0-based
    if (myrank <= maxrank0) then
      ng = ng0
      ig0 = myrank*ng0 ! 0-based index
      ig1 = ig0+ng0-1
    else !myrank > maxrank0
      ng = ng1
      ig0 = np0*ng0+(myrank-np0)*ng1 ! 0-based index
      ig1 = ig0+ng1-1
    endif
  else
    ! Error
    stop 'ERROR: illegal value of "np0"!'
  endif

  allocate(nzeros_row(ng))
  nzeros_row = 0
  do i = 1,nsparse1
    if (kgrow_sparse1(i)-1 >= ig0 .and. kgrow_sparse1(i)-1 <= ig1) then
      ind = kgrow_sparse1(i)-ig0 ! Fortran indexing
      nzeros_row(ind) = nzeros_row(ind)+1
    endif
  enddo
  !write(char_myrank,'(i4)') myrank

  call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,ngdof1,xvec1,ierr); CHECK_PETSC_ERROR(ierr)

  call VecDuplicate(xvec1,bvec1,ierr); CHECK_PETSC_ERROR(ierr)
  call VecDuplicate(xvec1,uvec1,ierr); CHECK_PETSC_ERROR(ierr)
  call VecDuplicate(xvec1,nzeror_gvec1,ierr); CHECK_PETSC_ERROR(ierr)
  call VecDuplicate(xvec1,nzeror_dvec1,ierr); CHECK_PETSC_ERROR(ierr)
  call VecDuplicate(xvec1,nzeror_ovec1,ierr); CHECK_PETSC_ERROR(ierr)
  call VecDuplicate(xvec1,iproc_gvec1,ierr); CHECK_PETSC_ERROR(ierr)
  call VecDuplicate(xvec1,interface_gvec1,ierr); CHECK_PETSC_ERROR(ierr)
  call VecDuplicate(xvec1,nself_gvec1,ierr); CHECK_PETSC_ERROR(ierr)
  call VecDuplicate(xvec1,ninterface_dvec1,ierr); CHECK_PETSC_ERROR(ierr)
  call VecDuplicate(xvec1,ninterface_ovec1,ierr); CHECK_PETSC_ERROR(ierr)

  if (SIMULATION_TYPE == 3) then
    ! Create backward xvector and associated b vec:
    call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,ngdof1,b_xvec1,ierr); CHECK_PETSC_ERROR(ierr)
    call VecDuplicate(b_xvec1,b_bvec1,ierr); CHECK_PETSC_ERROR(ierr)
  endif

  ! local vector
  call VecCreateSeq(PETSC_COMM_SELF,neq1,local_vec1,ierr)
  if (SIMULATION_TYPE == 3) then
    call VecCreateSeq(PETSC_COMM_SELF,neq1,b_local_vec1,ierr)
  endif

  ! objects needed for global vector scattering to local vector
  ! create local and global IS (index set) objects from the array of local and
  ! global indices
  call ISCreateGeneral(PETSC_COMM_WORLD,neq1,l2gdof1(1:),PETSC_COPY_VALUES,global_is,ierr); CHECK_PETSC_ERROR(ierr)
  call ISCreateStride(PETSC_COMM_SELF,neq1,0,1,local_is,ierr); CHECK_PETSC_ERROR(ierr)

  if (SIMULATION_TYPE == 3) then
    call ISCreateGeneral(PETSC_COMM_WORLD,neq1,l2gdof1(1:),PETSC_COPY_VALUES,b_global_is,ierr); CHECK_PETSC_ERROR(ierr)
    call ISCreateStride(PETSC_COMM_SELF,neq1,0,1,b_local_is,ierr); CHECK_PETSC_ERROR(ierr)
  endif

  ! create VecScatter object which is needed to scatter PETSc parallel vectors
  call VecScatterCreate(bvec1,global_is,local_vec1,local_is,vscat1,ierr); CHECK_PETSC_ERROR(ierr)

  if (SIMULATION_TYPE == 3) then
    call VecScatterCreate(b_bvec1, b_global_is, b_local_vec1, b_local_is, b_vscat1, ierr); CHECK_PETSC_ERROR(ierr)
  endif

  call ISDestroy(global_is,ierr) ! no longer necessary
  call ISDestroy(local_is,ierr)  ! no longer necessary

  if (SIMULATION_TYPE == 3) then
    call ISDestroy(b_global_is,ierr) ! no longer necessary
    call ISDestroy(b_local_is,ierr)  ! no longer necessary
  endif

  ! assign owner processor ID to each gdof (or row)
  allocate(ig_array1(ng),rproc_array1(ng))
  ig_array1 = (/ (i,i = ig0,ig1) /)

  !rproc=real(myrank)
  rproc_array1 = real(myrank)

  call VecSetValues(iproc_gvec1,ng,ig_array1,rproc_array1,INSERT_VALUES,ierr); CHECK_PETSC_ERROR(ierr)

  deallocate(ig_array1,rproc_array1)

  call VecAssemblyBegin(iproc_gvec1,ierr); CHECK_PETSC_ERROR(ierr)
  call VecAssemblyEnd(iproc_gvec1,ierr); CHECK_PETSC_ERROR(ierr)

  call VecMin(iproc_gvec1,PETSC_NULL_INTEGER,pmin,ierr)
  call VecMax(iproc_gvec1,PETSC_NULL_INTEGER,pmax,ierr)

  if (myrank == 0) print *,'PETSc solver: iproc range',pmin,pmax
  call synchronize_all()

  ! copy solution to local array
  allocate(iproc_array1(neq1),rproc_array1(neq1))
  call scatter_globalvec1(iproc_gvec1, rproc_array1)

  iproc_array1 = int(rproc_array1)

  if (myrank == 0) print *,'PETSc solver: vector3 iproc',minval(iproc_array1),maxval(iproc_array1)
  call synchronize_all()

  !!TODO: use local scatter
  !call VecScatterCreateToAll(iproc_gvec1,pscat1,iproc_garray,ierr);
  !call VecScatterBegin(pscat1,iproc_gvec1,iproc_garray,INSERT_VALUES,SCATTER_FORWARD,ierr);
  !call VecScatterEnd(pscat1,iproc_gvec1,iproc_garray,INSERT_VALUES,SCATTER_FORWARD,ierr);
  !call VecScatterDestroy(pscat1);

  ! assign interface ID to each gdofs
  rval = 1.0
  ! inner core
  do i = 1,num_interfaces_inner_core1
    nibool = nibool_interfaces_inner_core1(i)
    allocate(ibool_interface(nibool))
    ibool_interface = ibool_interfaces_inner_core1(1:nibool,i)
    do i_bool = 1,nibool
      do i_ndof = 1,NNDOF
        igdof = ggdof_ic1(i_ndof,ibool_interface(i_bool))-1
        if (igdof >= 0) call VecSetValues(interface_gvec1,1,igdof,rval,INSERT_VALUES,ierr)
      enddo
    enddo
    deallocate(ibool_interface)
  enddo

  ! outer core
  do i = 1,num_interfaces_outer_core1
    nibool = nibool_interfaces_outer_core1(i)
    allocate(ibool_interface(nibool))
    ibool_interface = ibool_interfaces_outer_core1(1:nibool,i)
    do i_bool = 1,nibool
      do i_ndof = 1,NNDOF
        igdof = ggdof_oc1(i_ndof,ibool_interface(i_bool))-1
        if (igdof >= 0) call VecSetValues(interface_gvec1,1,igdof,rval,INSERT_VALUES,ierr)
      enddo
    enddo
    deallocate(ibool_interface)
  enddo

  ! crust mantle
  do i = 1,num_interfaces_crust_mantle1
    nibool = nibool_interfaces_crust_mantle1(i)
    allocate(ibool_interface(nibool))
    ibool_interface = ibool_interfaces_crust_mantle1(1:nibool,i)
    do i_bool = 1,nibool
      do i_ndof = 1,NNDOF
        igdof = ggdof_cm1(i_ndof,ibool_interface(i_bool))-1
        if (igdof >= 0) call VecSetValues(interface_gvec1,1,igdof,rval,INSERT_VALUES,ierr)
      enddo
    enddo
    deallocate(ibool_interface)
  enddo

  ! transition infinite
  if (ADD_TRINF) then
    do i = 1,num_interfaces_trinfinite1
      nibool = nibool_interfaces_trinfinite1(i)
      allocate(ibool_interface(nibool))
      ibool_interface = ibool_interfaces_trinfinite1(1:nibool,i)
      do i_bool = 1,nibool
        do i_ndof = 1,NNDOF
          igdof = ggdof_trinf1(i_ndof,ibool_interface(i_bool))-1
          if (igdof >= 0) call VecSetValues(interface_gvec1,1,igdof,rval,INSERT_VALUES,ierr)
        enddo
      enddo
      deallocate(ibool_interface)
    enddo
  endif

  ! infinite
  do i = 1,num_interfaces_infinite1
    nibool = nibool_interfaces_infinite1(i)
    allocate(ibool_interface(nibool))
    ibool_interface = ibool_interfaces_infinite1(1:nibool,i)
    do i_bool = 1,nibool
      do i_ndof = 1,NNDOF
        igdof = ggdof_inf1(i_ndof,ibool_interface(i_bool))-1
        if (igdof >= 0) call VecSetValues(interface_gvec1,1,igdof,rval,INSERT_VALUES,ierr)
      enddo
    enddo
    deallocate(ibool_interface)
  enddo

  call VecAssemblyBegin(interface_gvec1,ierr); CHECK_PETSC_ERROR(ierr)
  call VecAssemblyEnd(interface_gvec1,ierr); CHECK_PETSC_ERROR(ierr)

  ! copy solution to local array
  allocate(isg_interface(neq1),rg_interface(neq1))
  call scatter_globalvec1(interface_gvec1,rg_interface)
  isg_interface = int(rg_interface)

  ! estimate correction for the number of nonzero entries in the diagonal and nondiagonal portion
  ! self interface
  !rval=-1.0
  !call VecSet(nself_gvec1,rval,ierr) ! subtract self

  rval = 1.0
  do i = 1,neq1
    if (isg_interface(i) == 1) then
      call VecSetValues(nself_gvec1,1,l2gdof1(i),rval,ADD_VALUES,ierr);
    endif
  enddo
  call VecAssemblyBegin(nself_gvec1,ierr); CHECK_PETSC_ERROR(ierr)
  call VecAssemblyEnd(nself_gvec1,ierr); CHECK_PETSC_ERROR(ierr)
  call VecGetLocalSize(nself_gvec1,n,ierr)

  allocate(rnself_lgarray1(neq1))
  call scatter_globalvec1(nself_gvec1, rnself_lgarray1)
  call VecGetArrayF90(nself_gvec1,rnself_array1,ierr)

  allocate(nself_array1(n))
  nself_array1 = int(rnself_array1(1:n))

  where(nself_array1 > 0) nself_array1 = nself_array1-1 ! subtract self

  call VecRestoreArrayF90(nself_gvec1,rnself_array1,ierr)
  call VecDestroy(nself_gvec1,ierr)

  if (myrank == 0) print *,'PETSc solver: maximum value of nself:',maxval(nself_array1)

  ! factor for maximum number of interfaces for each nondiagonal entry of the stiffness matrix
  ! the factor below is valid ONLY for rectangular partitioning of the global model
  max_ni = 8.0
  fac_ni = 0.0

  ! first element
  igr0 = kgrow_sparse1(1)-1
  igc0 = kgcol_sparse1(1)-1
  ir0 = krow_sparse1(1)
  ic0 = kcol_sparse1(1)
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
    igr = kgrow_sparse1(i)-1
    igc = kgcol_sparse1(i)-1
    ir = krow_sparse1(i)
    ic = kcol_sparse1(i)
    if (l2gdof1(ir) /= igr .or. l2gdof1(ic) /= igc) then
      print *,'Error: strange:',l2gdof1(ir),igr,l2gdof1(ic),igc
      stop
    endif
    if (igr /= igr0) then
      ! new row starts
      ! set values computed so far
      rnd = real(nd)
      rnoffd = real(noffd)
      call VecSetValues(nzeror_dvec1,1,igr0,rnd,ADD_VALUES,ierr); CHECK_PETSC_ERROR(ierr)
      call VecSetValues(nzeror_ovec1,1,igr0,rnoffd,ADD_VALUES,ierr); CHECK_PETSC_ERROR(ierr)

      call VecSetValues(ninterface_dvec1,1,igr0,rnid,ADD_VALUES,ierr); CHECK_PETSC_ERROR(ierr)
      call VecSetValues(ninterface_ovec1,1,igr0,rnioffd,ADD_VALUES,ierr); CHECK_PETSC_ERROR(ierr)

      ! reset
      nd = 0; noffd = 0
      rnid = 0.; rnioffd = 0.
      igr0 = igr !kgrow_sparse1(i)-1
      igc0 = igc !kgcol_sparse1(i)-1
      ir0 = ir !krow_sparse1(i)
      ic0 = ic !kcol_sparse1(i)

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
          rnid = rnid+(1.0/fac_ni)
        endif
      else
        noffd = noffd+1
        if (igr /= igc .and. rnself_lgarray1(ir) > 0.0 .and. rnself_lgarray1(ic) > 0.0) then
          fac_ni = min(max_ni,min(rnself_lgarray1(ir),rnself_lgarray1(ic)))
          rnioffd = rnioffd+(1.0/fac_ni)
        endif
      endif
    endif
    if (i == nsparse1) then
      ! for last
      rnd = real(nd)
      rnoffd = real(noffd)
      call VecSetValues(nzeror_dvec1,1,igr0,rnd,ADD_VALUES,ierr); CHECK_PETSC_ERROR(ierr)
      call VecSetValues(nzeror_ovec1,1,igr0,rnoffd,ADD_VALUES,ierr); CHECK_PETSC_ERROR(ierr)

      call VecSetValues(ninterface_dvec1,1,igr0,rnid,ADD_VALUES,ierr); CHECK_PETSC_ERROR(ierr)
      call VecSetValues(ninterface_ovec1,1,igr0,rnioffd,ADD_VALUES,ierr); CHECK_PETSC_ERROR(ierr)
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
  ! diagonal matrix
  call VecGetLocalSize(nzeror_dvec1,n,ierr)
  call VecGetArrayF90(nzeror_dvec1,nzeror_darray1,ierr)

  allocate(nnzero_diag1(n))
  nnzero_diag1(:) = int(nzeror_darray1(1:n))
  nnzero_diag1(:) = nnzero_diag1(:) - nself_array1(:)

  !debug
  if (myrank == 0) print *,'PETSc solver: n: ',n,minval(nzeror_darray1),maxval(nzeror_darray1), &
                                                 minval(nnzero_diag1),maxval(nnzero_diag1)
  call synchronize_all()

  call VecRestoreArrayF90(nzeror_dvec1,nzeror_darray1,ierr)
  call VecDestroy(nzeror_dvec1,ierr)

  ! off-diagonal matrix
  call VecGetArrayF90(nzeror_ovec1,nzeror_oarray1,ierr)

  allocate(nnzero_offdiag1(n))
  nnzero_offdiag1(:) = int(nzeror_oarray1(1:n))

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

  where(ninterface_darray1 > 0) ninterface_darray1 = ninterface_darray1 - 4
  where(ninterface_darray1 < 0) ninterface_darray1 = 0

  call VecGetArrayF90(ninterface_ovec1,rninterface_oarray1,ierr)

  !where(rninterface_oarray1>0.0 .and. rninterface_oarray1 < 1.0)rninterface_oarray1=1.0

  allocate(ninterface_oarray1(n))
  ninterface_oarray1 = int(rninterface_oarray1(1:n))

  call VecRestoreArrayF90(ninterface_ovec1,rninterface_oarray1,ierr)
  call VecDestroy(ninterface_ovec1,ierr)

  where(ninterface_oarray1 > 0) ninterface_oarray1 = ninterface_oarray1 - 8
  where(ninterface_oarray1 < 0) ninterface_oarray1 = 0

  nnzero_diag1(:) = nnzero_diag1(:) - ninterface_darray1(:)
  nnzero_offdiag1(:) = nnzero_offdiag1(:) - ninterface_oarray1(:)

  !debug
  if (myrank == 0) print *,'PETSc solver: non-zero diag   :',minval(nnzero_diag1),maxval(nnzero_diag1)
  if (myrank == 0) print *,'PETSc solver: non-zero offdiag:',minval(nnzero_offdiag1),maxval(nnzero_offdiag1)
  call synchronize_all()

  rval = 1.0
  do i = 1,nsparse1
    igdof = kgrow_sparse1(i)-1 ! Fortran index
    call VecSetValues(nzeror_gvec1,1,igdof,rval,ADD_VALUES,ierr); CHECK_PETSC_ERROR(ierr)
  enddo
  call VecAssemblyBegin(nzeror_gvec1,ierr); CHECK_PETSC_ERROR(ierr)
  call VecAssemblyEnd(nzeror_gvec1,ierr); CHECK_PETSC_ERROR(ierr)

  call VecGetLocalSize(nzeror_gvec1,n,ierr); CHECK_PETSC_ERROR(ierr)

  !debug
  if (myrank == 0) print *,'PETSc solver: size of vector:',ng,n,minval(kgrow_sparse1),ig0
  call synchronize_all()

  ! free temporary array
  deallocate(kgrow_sparse1,kgcol_sparse1)

  ! non-zero array for diagonal/off-diagonal matrix?
  call VecGetArrayF90(nzeror_gvec1,nzeror_array1,ierr); CHECK_PETSC_ERROR(ierr)

  allocate(inzeror_array1(n))
  inzeror_array1(:) = int(nzeror_array1(1:n))

  call VecRestoreArrayF90(nzeror_gvec1,nzeror_array1,ierr); CHECK_PETSC_ERROR(ierr)
  call VecDestroy(nzeror_gvec1,ierr); CHECK_PETSC_ERROR(ierr)

  inzeros_max = maxvec(inzeror_array1)
  inzeros_min = minvec(inzeror_array1)

  ! Create the stiffness matrix (same for forward/adjoint simulations)
  call MatCreate(PETSC_COMM_WORLD,Amat1,ierr); CHECK_PETSC_ERROR(ierr)
  call MatSetType(Amat1,MATMPIAIJ,ierr); CHECK_PETSC_ERROR(ierr)
  call MatSetSizes(Amat1,PETSC_DECIDE,PETSC_DECIDE,ngdof1,ngdof1,ierr); CHECK_PETSC_ERROR(ierr)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    number of global degrees of freedom   : ',ngdof1
    write(IMAIN,*) '    number of non-zero row entries min/max: ',nzeros_min,'/',nzeros_max
    write(IMAIN,*) '    number of inzeror entries min/max     : ',inzeros_min,'/',inzeros_max
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! preallocation
  !TODO: please re-asses this preallocation for different meshes and/or correct the diagonal/off-diagonal non-zero row entries
  !
  ! way 1: explicitly specifies diagonal and off-diagonal local submatrices by setting array inzeror_array1,
  !        containing the number of non-zeros per row.
  !
  !        this will lead to a lot of re-allocations when setting the region values by MatSetValues() in petsc_set_matrix1(),
  !        as the non-zero estimates seem to be off.
  !
  !call MatMPIAIJSetPreallocation(Amat1,nzeros_max,inzeror_array1,nzeros_max,inzeror_array1,ierr); CHECK_PETSC_ERROR(ierr)

  ! way 2: only specify the nonzero structure by nzeros_max, and let Petsc decide about array structuring.
  !
  !        this seems to lead to a much faster petsc_set_matrix1() routine without the re-allocations.
  !        however, the diagonal and in particular the off-diagonal estimate with nzeros_max might be still off.
  call MatMPIAIJSetPreallocation(Amat1,nzeros_max,PETSC_NULL_INTEGER, &
                                 nzeros_max,PETSC_NULL_INTEGER,ierr); CHECK_PETSC_ERROR(ierr)

  call MatSetFromOptions(Amat1,ierr); CHECK_PETSC_ERROR(ierr)
  call MatGetOwnershipRange(Amat1,istart,iend,ierr); CHECK_PETSC_ERROR(ierr)

  ! check
  if (istart /= ig0 .or. iend-1 /= ig1) then
    print *,'ERROR: PETSc solver ownership range mismatch!'
    print *,'       ownership range:',myrank,istart,ig0,iend-1,ig1,nzeros_row(1)
    stop
  endif
  call synchronize_all()

  ! free temporary array
  deallocate(nzeros_row)

  ! memory usage (returns bytes)
  call PetscMemoryGetCurrentUsage(bytes, ierr); CHECK_PETSC_ERROR(ierr)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    current PETSc memory usage  = ',sngl(bytes / 1024.d0 / 1024.d0),'MB'
    write(IMAIN,*) '                                = ',sngl(bytes / 1024.d0 / 1024.d0 / 1024.d0),'GB'
    call flush_IMAIN()
  endif

  ! Create forward solver
  !solver_type1 = CG
  !call create_linear_solver(solver_type1, ksp1, Amat1, pc1, Fmat1)

  !call KSPSetTolerances(ksp1,KSP_RTOL1,KSP_ATOL1,KSP_DTOL1,KSP_MAXITER1,ierr); CHECK_PETSC_ERROR(ierr)

  !  Set runtime options, e.g.,
  !    -ksp_type < type> -pc_type < type> -ksp_monitor -ksp_KSP_RTOL < KSP_RTOL>
  !  These options will override those specified above as long as
  !  KSPSetFromOptions() is called _after_ any other customization
  !  routines.
  !call KSPSetFromOptions(ksp1,ierr)

  ! Create adjoint solver:
  !if (SIMULATION_TYPE == 3) then
  !  call create_linear_solver(solver_type1, b_ksp1, Amat1, pc1, Fmat1)
  !  call KSPSetTolerances(b_ksp1,KSP_RTOL1,KSP_ATOL1,KSP_DTOL1,KSP_MAXITER1,ierr); CHECK_PETSC_ERROR(ierr)
  !  call KSPSetFromOptions(b_ksp1,ierr)
  !  if (myrank == 0) then
  !    print *,'PETSc solver: Created adjoint linear KSP solver...'
  !  endif
  !endif

  !-------------------------------------------------------------------------------
  ! Create the linear solver and set various options
  !-------------------------------------------------------------------------------
  ! define solver type
  ! COMMAND: define from the command
  ! SUPERLU: SuperLU solver
  ! MUMPS: MUMPS solver
  !solver_type1 = CG
  ! Create linear solver context

  call KSPCreate(PETSC_COMM_WORLD,ksp1,ierr); CHECK_PETSC_ERROR(ierr)
  call KSPSetOperators(ksp1,Amat1,Amat1,ierr); CHECK_PETSC_ERROR(ierr) ! version >= 3.5

  if (SIMULATION_TYPE == 3) then
    call KSPSetInitialGuessNonzero(ksp1,PETSC_FALSE,ierr); CHECK_PETSC_ERROR(ierr)
  else
    call KSPSetInitialGuessNonzero(ksp1,PETSC_TRUE,ierr); CHECK_PETSC_ERROR(ierr)
  endif

  call KSPSetDiagonalScale(ksp1,PETSC_TRUE,ierr); CHECK_PETSC_ERROR(ierr)
  call KSPSetReusePreconditioner(ksp1,PETSC_TRUE,ierr); CHECK_PETSC_ERROR(ierr)
  call KSPSetType(ksp1,KSPCG,ierr); CHECK_PETSC_ERROR(ierr)
  call KSPGetPC(ksp1,pc1,ierr); CHECK_PETSC_ERROR(ierr)
  call PCFactorSetShiftType(pc1,MAT_SHIFT_POSITIVE_DEFINITE,ierr); CHECK_PETSC_ERROR(ierr)
  call KSPSetTolerances(ksp1,KSP_RTOL1,KSP_ATOL1,KSP_DTOL1,KSP_MAXITER1,ierr); CHECK_PETSC_ERROR(ierr)
  call KSPSetFromOptions(ksp1,ierr); CHECK_PETSC_ERROR(ierr)

  ! BACKWARD SOLVER
  if (SIMULATION_TYPE == 3) then
    call KSPCreate(PETSC_COMM_WORLD,b_ksp1,ierr); CHECK_PETSC_ERROR(ierr)
    call KSPSetOperators(b_ksp1,Amat1,Amat1,ierr); CHECK_PETSC_ERROR(ierr) ! version >= 3.5
    call KSPSetInitialGuessNonzero(b_ksp1,PETSC_FALSE,ierr); CHECK_PETSC_ERROR(ierr)
    call KSPSetDiagonalScale(b_ksp1,PETSC_TRUE,ierr); CHECK_PETSC_ERROR(ierr)
    call KSPSetReusePreconditioner(b_ksp1,PETSC_TRUE,ierr); CHECK_PETSC_ERROR(ierr)
    call KSPSetType(b_ksp1,KSPCG,ierr); CHECK_PETSC_ERROR(ierr)
    call KSPGetPC(b_ksp1,b_pc1,ierr); CHECK_PETSC_ERROR(ierr)
    call PCFactorSetShiftType(b_pc1,MAT_SHIFT_POSITIVE_DEFINITE,ierr); CHECK_PETSC_ERROR(ierr)
    call KSPSetTolerances(b_ksp1,KSP_RTOL1,KSP_ATOL1,KSP_DTOL1,KSP_MAXITER1,ierr); CHECK_PETSC_ERROR(ierr)
    call KSPSetFromOptions(b_ksp1,ierr); CHECK_PETSC_ERROR(ierr)
  endif

  !debug
  if (myrank == 0) print *,'PETSc solver: ---------- Finished PETSC initialisation ---------- '
  if (myrank == 0) print *

  ! synchronize all processes
  call synchronize_all()

  ! user output
  if (myrank == 0) then
    tCPU = wtime() - tstart
    write(IMAIN,*) '    elapsed time for PETSc solver initialization: ',sngl(tCPU),'(s)'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

#else
  ! no PETSc compilation support
  ! compilation without PETSc support
  if (myrank == 0) then
    print *, "Error: PETSc solver enabled without PETSc Support."
    print *, "       To enable PETSc support, reconfigure with --with-petsc flag."
  endif
  ! safety stop
  call exit_MPI(myrank,"Error PETSc solver: petsc_initialize1 called without compilation support")

#endif

  end subroutine petsc_initialize1

!
!===============================================================================
!

#ifdef USE_PETSC

! not used so far...

!  subroutine create_linear_solver(stype, l_ksp, l_Amat, l_pc, l_fmat)
!
!  ! Create the linear solver and set various options
!  ! stype: Solver type - options available are
!  !   COMMAND   define from the command
!  !   SUPERLU   SuperLU solver
!  !   MUMPS     MUMPS solver
!  ! l_ksp: the local KSP (ksp1 or b_ksp1 etc)
!  ! l_Amat: local A matrix - I think always Amat1
!  ! l_pc:   local preconditioner e.g pc1
!
!  use specfem_par, only: myrank, SIMULATION_TYPE
!
!  implicit none
!
!  PetscInt :: stype
!  type(tKSP) :: l_ksp
!  type(tMat) :: l_Amat, l_fmat
!  type(tPC) :: l_pc
!
!  ! Create linear solver context
!  call KSPCreate(PETSC_COMM_WORLD,l_ksp,ierr)
!  ! Set operators. Here the matrix that defines the linear system
!  ! also serves as the preconditioning matrix.
!  !call KSPSetOperators(ksp1,Amat1,Amat1,SAME_PRECONDITIONER,ierr) ! version < 3.5
!  call KSPSetOperators(l_ksp,l_Amat,l_Amat,ierr) ! version >= 3.5
!
!  call KSPSetInitialGuessNonzero(l_ksp,PETSC_TRUE,ierr); CHECK_PETSC_ERROR(ierr)
!  !since the euqutions are nondimensionalized, the scaling is unnecessary?
!  call KSPSetDiagonalScale(l_ksp,PETSC_TRUE,ierr); CHECK_PETSC_ERROR(ierr)
!  call KSPSetReusePreconditioner(l_ksp,PETSC_TRUE,ierr)
!
!  if (stype == COMMAND) then
!    if (myrank == 0) print *,'Solver type: provided via command'
!  else if (stype == CG) then
!    ! CONJUGATE GRADIENT
!    if (myrank == 0) print *,'Solver type: CG'
!    call KSPSetType(l_ksp,KSPCG,ierr); CHECK_PETSC_ERROR(ierr)
!    ! Fetch preconditioner
!    call KSPGetPC(l_ksp,l_pc,ierr); CHECK_PETSC_ERROR(ierr)
!    call PCFactorSetShiftType(l_pc,MAT_SHIFT_POSITIVE_DEFINITE,ierr); CHECK_PETSC_ERROR(ierr)
!  else if (stype == SUPERLU) then
!    ! SUPER LU
!    if (myrank == 0) print *,'Solver type: SUPERLU'
!    if (SIMULATION_TYPE == 3) then
!      stop ' ERROR: SUPERLU not implemented for adjoint sims yet.'
!    endif
!    flg_ilu    = PETSC_FALSE;
!    flg_lu     = PETSC_FALSE;
!    ! version < 3.8.0
!    !call PetscOptionsGetBool(PETSC_NULL_CHARACTER,"-use_superlu_lu",flg_lu,flg,ierr);
!    call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
!          "-use_superlu_lu",flg_lu,flg,ierr); CHECK_PETSC_ERROR(ierr)
!    !call PetscOptionsGetBool(PETSC_NULL_CHARACTER,"-use_superlu_ilu",flg_ilu,flg,ierr);
!    if (flg_lu .or. flg_ilu) then
!      call KSPSetType(l_ksp,KSPPREONLY,ierr); CHECK_PETSC_ERROR(ierr)
!      call KSPGetPC(l_ksp,l_pc,ierr); CHECK_PETSC_ERROR(ierr)
!      if (flg_lu) then
!        call PCSetType(l_pc,PCLU,ierr); CHECK_PETSC_ERROR(ierr)
!      else if (flg_ilu) then
!        call PCSetType(l_pc,PCILU,ierr); CHECK_PETSC_ERROR(ierr)
!      endif
!      call PCFactorSetShiftType(l_pc,MAT_SHIFT_POSITIVE_DEFINITE,ierr); CHECK_PETSC_ERROR(ierr)
!      ! version < 3.9
!      !call PCFactorSetMatSolverPackage(l_pc,MATSOLVERSUPERLU,ierr)
!      call PCFactorSetMatSolverType(l_pc,MATSOLVERSUPERLU,ierr); CHECK_PETSC_ERROR(ierr)
!      ! version < 3.9
!      !call PCFactorSetUpMatSolverPackage(l_pc,ierr); ! call MatGetFactor() to create F
!      call PCFactorSetUpMatSolverType(l_pc,ierr); CHECK_PETSC_ERROR(ierr)
!      call PCFactorGetMatrix(l_pc,l_fmat,ierr); CHECK_PETSC_ERROR(ierr)
!      !call MatSuperluSetILUDropTol(l_fmat,1.e-8,ierr); CHECK_PETSC_ERROR(ierr)
!    endif
!  else if (stype == MUMPS) then
!    if (myrank == 0) print *,'Solver type: MUMPS'
!
!    stop 'ERROR - WE commented out MUMPS stuff due to syntax error'
!
!    flg_lu    = PETSC_FALSE;
!    flg_ch    = PETSC_FALSE;
!    ! version < 3.8.0
!    ! call PetscOptionsGetBool(PETSC_NULL_CHARACTER,"-use_mumps_ch",flg_ch,flg,ierr);
!     call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
!          "-use_mumps_ch",flg_ch,flg,ierr);
!    if (flg_lu .or. flg_ch) then
!      call KSPSetType(l_ksp,KSPPREONLY,ierr); CHECK_PETSC_ERROR(ierr)
!      call KSPGetPC(l_ksp,l_pc,ierr); CHECK_PETSC_ERROR(ierr)
!      if (flg_lu) then
!        call PCSetType(l_pc,PCLU,ierr); CHECK_PETSC_ERROR(ierr)
!      else if (flg_ch) then
!        ! set MUMPS id%SYM=1
!        call MatSetOption(l_Amat,MAT_SPD,PETSC_TRUE,ierr); CHECK_PETSC_ERROR(ierr)
!        call PCSetType(l_pc,PCCHOLESKY,ierr);
!      endif
!      call PCFactorSetShiftType(l_pc,MAT_SHIFT_POSITIVE_DEFINITE,ierr); CHECK_PETSC_ERROR(ierr)
!      ! version < 3.9
!      !call PCFactorSetMatSolverPackage(l_pc,MATSOLVERMUMPS,ierr);
!      !call PCFactorSetUpMatSolverPackage(l_pc,ierr); ! call MatGetFactor() to create F
!      !call PCFactorSetMatSolverType(l_pc,MATSOLVERMUMPS,ierr);
!      !call PCFactorSetUpMatSolverType(l_pc,ierr); ! call MatGetFactor() to create F
!      !call PCFactorGetMatrix(l_pc,l_fmat,ierr);
!      icntl = 7; ival = 2;
!      !call MatMumpsSetIcntl(l_fmat,icntl,ival,ierr);
!      icntl = 1; val = 0.0;
!      !call MatMumpsSetCntl(l_fmat,icntl,val,ierr);
!    endif
!  endif
!
!  end subroutine create_linear_solver

#endif

!
!===============================================================================
!

  subroutine petsc_set_matrix1()

#ifdef USE_PETSC

  use constants, only: CUSTOM_REAL,IFLAG_IN_FICTITIOUS_CUBE,NGLLCUBE_INF

  use specfem_par, only: NSPEC_INNER_CORE, &
    NSPEC_OUTER_CORE,NSPEC_CRUST_MANTLE,NSPEC_TRINFINITE,NSPEC_INFINITE

  use specfem_par_innercore, only: idoubling_inner_core

  use specfem_par_full_gravity, only: NEDOF,NEDOF1, &
    ggdof_ic1,storekmat_inner_core1,inode_elmt_ic1, &
    ggdof_oc1,storekmat_outer_core1,inode_elmt_oc1, &
    ggdof_cm1,storekmat_crust_mantle1,inode_elmt_cm1, &
    ggdof_trinf1,storekmat_trinfinite1,inode_elmt_trinf1, &
    ggdof_inf1,storekmat_infinite1,inode_elmt_inf1

  implicit none
  integer :: i,i_elmt,j !,ncount
  integer :: ggdof_elmt(NEDOF1) !,idof(NEDOF1),igdof(NEDOF1)

  PetscInt:: istart,iend,ndiag,noffdiag

  !debugging
  logical, parameter :: DEBUG_FILE_OUTPUT = .false.
  integer :: ncols
  integer,dimension(:),allocatable :: cols

  character(len=10) :: char_myrank
  character(len=60) :: outf_name

  ! types required by MatSetValues:
  !   MatSetValues(Mat mat, PetscInt m, const PetscInt idxm[], PetscInt n, \
  !                const PetscInt idxn[], const PetscScalar v[], InsertMode addv)
  PetscScalar :: v

  ! matrix info
  double precision :: info(MAT_INFO_SIZE)
  double precision :: mallocsval
  PetscLogDouble :: bytes

  ! timing
  double precision :: tstart,tCPU
  double precision, external :: wtime

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    setting up Level-1 solver matrix'
    call flush_IMAIN()
  endif

  ! timing
  tstart = wtime()

  call MatZeroEntries(Amat1,ierr); CHECK_PETSC_ERROR(ierr)

  ! to avoid errors:
  !   note that the preallocation of Amat1 might be slightly off, in that the number of non-zero elements
  !   for the diagonal and off-diagonal matrices could be estimated wrongly.
  !
  !   assigning the matrix entry values below with MatSetValues() might lead to errors like:
  !       [4]PETSC ERROR: Argument out of range
  !       [4]PETSC ERROR: New nonzero at (11346,11355) caused a malloc
  !       Use MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE) to turn off this check
  !       [4]PETSC ERROR: See https://petsc.org/release/faq/ for trouble shooting.
  !       ..
  !   and abort execution.
  !
  ! this is to turn off new malloc check error messages, when a new malloc is required by MatSetValues()
  !call MatSetOption(Amat1, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr); CHECK_PETSC_ERROR(ierr)

  ! note: although the execution works by turning this check off, the matrix value assignment in particular
  !       for the crust/mantle region below takes very long. this might be due to excessive malloc's required.
  !
  ! TODO: it would be great to fix the preallocation of Amat1 and specify more exact non-zero matrix entries
  !       for anybody who knows how to do this :)

  ! user info
  if (myrank == 0) then
    write(IMAIN,*) '    setup inner core values'
    call flush_IMAIN()
  endif

  ! inner core
  do i_elmt = 1,NSPEC_INNER_CORE
    if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE) cycle

    ggdof_elmt(:) = reshape(ggdof_ic1(:,inode_elmt_ic1(:,i_elmt)),(/NEDOF1/))
    ggdof_elmt(:) = ggdof_elmt(:)-1 ! petsc index starts from 0
    !ncount = 0; idof(:) = -1; igdof(:) = -1
    do i = 1,NEDOF1     ! NGLLCUBE_INF * NNDOF
      do j = 1,NEDOF1
        if (ggdof_elmt(i) >= 0 .and. ggdof_elmt(j) >= 0 .and. &
            storekmat_inner_core1(i,j,i_elmt) /= 0.0_CUSTOM_REAL) then
          !call MatSetValues(Amat1,1,ggdof_elmt(i),1,ggdof_elmt(j), &
          !                  storekmat_inner_core1(i,j,i_elmt),ADD_VALUES,ierr); CHECK_PETSC_ERROR(ierr)
          ! petsc types
          v = storekmat_inner_core1(i,j,i_elmt)
          call MatSetValues(Amat1,1,ggdof_elmt(i),1,ggdof_elmt(j),v,ADD_VALUES,ierr); CHECK_PETSC_ERROR(ierr)
        endif
      enddo
    enddo
  enddo
  call synchronize_all()
  deallocate(storekmat_inner_core1)

  ! user output
  if (myrank == 0) then
    tCPU = wtime() - tstart
    write(IMAIN,*) '    elapsed time: ',sngl(tCPU),'(s)'
    write(IMAIN,*) '    setup outer core values'
    call flush_IMAIN()
  endif

  ! outer core
  do i_elmt = 1,NSPEC_OUTER_CORE
    ggdof_elmt(:) = reshape(ggdof_oc1(:,inode_elmt_oc1(:,i_elmt)),(/NEDOF1/))
    ggdof_elmt(:) = ggdof_elmt(:)-1 ! petsc index starts from 0
    !ncount = 0; idof(:) = -1; igdof(:) = -1
    do i = 1,NEDOF1
      do j = 1,NEDOF1
        if (ggdof_elmt(i) >= 0 .and. ggdof_elmt(j) >= 0 .and. &
            storekmat_outer_core1(i,j,i_elmt) /= 0.0_CUSTOM_REAL) then
          !call MatSetValues(Amat1,1,ggdof_elmt(i),1,ggdof_elmt(j), &
          !                  storekmat_outer_core1(i,j,i_elmt),ADD_VALUES,ierr); CHECK_PETSC_ERROR(ierr)
          ! petsc types
          v = storekmat_outer_core1(i,j,i_elmt)
          call MatSetValues(Amat1,1,ggdof_elmt(i),1,ggdof_elmt(j),v,ADD_VALUES,ierr); CHECK_PETSC_ERROR(ierr)
        endif
      enddo
    enddo
  enddo
  call synchronize_all()
  deallocate(storekmat_outer_core1)

  ! user output
  if (myrank == 0) then
    tCPU = wtime() - tstart
    write(IMAIN,*) '    elapsed time: ',sngl(tCPU),'(s)'
    write(IMAIN,*) '    setup crust/mantle values'
    call flush_IMAIN()
  endif

  ! crust mantle
  do i_elmt = 1,NSPEC_CRUST_MANTLE
    ggdof_elmt(:) = reshape(ggdof_cm1(:,inode_elmt_cm1(:,i_elmt)),(/NEDOF1/))
    ggdof_elmt(:) = ggdof_elmt(:)-1 ! petsc index starts from 0
    !ncount = 0; idof(:) = -1; igdof(:) = -1
    do i = 1,NEDOF1
      do j = 1,NEDOF1
        if (ggdof_elmt(i) >= 0 .and. ggdof_elmt(j) >= 0 .and. &
            storekmat_crust_mantle1(i,j,i_elmt) /= 0.0_CUSTOM_REAL) then
          !call MatSetValues(Amat1,1,ggdof_elmt(i),1,ggdof_elmt(j), &
          !                  storekmat_crust_mantle1(i,j,i_elmt),ADD_VALUES,ierr); CHECK_PETSC_ERROR(ierr)
          ! petsc types
          v = storekmat_crust_mantle1(i,j,i_elmt)
          call MatSetValues(Amat1,1,ggdof_elmt(i),1,ggdof_elmt(j),v,ADD_VALUES,ierr); CHECK_PETSC_ERROR(ierr)
        endif
      enddo
    enddo
  enddo
  call synchronize_all()
  deallocate(storekmat_crust_mantle1)

  ! user output
  if (myrank == 0) then
    tCPU = wtime() - tstart
    write(IMAIN,*) '    elapsed time: ',sngl(tCPU),'(s)'
    call flush_IMAIN()
  endif

  ! trinfinite
  if (NSPEC_INFINITE > 0) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '    setup trinfinite values'
      call flush_IMAIN()
    endif

    do i_elmt = 1,NSPEC_TRINFINITE
      ggdof_elmt(:) = reshape(ggdof_trinf1(:,inode_elmt_trinf1(:,i_elmt)),(/NEDOF1/))
      ggdof_elmt(:) = ggdof_elmt(:)-1 ! petsc index starts from 0
      !ncount = 0; idof(:) = -1; igdof(:) = -1
      do i = 1,NEDOF1
        do j = 1,NEDOF1
          if (ggdof_elmt(i) >= 0 .and. ggdof_elmt(j) >= 0 .and. &
              storekmat_trinfinite1(i,j,i_elmt) /= 0.0_CUSTOM_REAL) then
            !call MatSetValues(Amat1,1,ggdof_elmt(i),1,ggdof_elmt(j), &
            !                  storekmat_trinfinite1(i,j,i_elmt),ADD_VALUES,ierr); CHECK_PETSC_ERROR(ierr)
            ! petsc types
            v = storekmat_trinfinite1(i,j,i_elmt)
            call MatSetValues(Amat1,1,ggdof_elmt(i),1,ggdof_elmt(j),v,ADD_VALUES,ierr); CHECK_PETSC_ERROR(ierr)
          endif
        enddo
      enddo
    enddo
    ! user output
    if (myrank == 0) then
      tCPU = wtime() - tstart
      write(IMAIN,*) '    elapsed time: ',sngl(tCPU),'(s)'
    endif
  endif
  call synchronize_all()
  deallocate(storekmat_trinfinite1)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    setup infinite values'
    call flush_IMAIN()
  endif

  ! infinite
  do i_elmt = 1,NSPEC_INFINITE
    ggdof_elmt(:) = reshape(ggdof_inf1(:,inode_elmt_inf1(:,i_elmt)),(/NEDOF1/))
    ggdof_elmt(:) = ggdof_elmt(:)-1 ! petsc index starts from 0
    !ncount = 0; idof(:) = -1; igdof(:) = -1
    do i = 1,NEDOF1
      do j = 1,NEDOF1
        if (ggdof_elmt(i) >= 0 .and. ggdof_elmt(j) >= 0 .and. &
            storekmat_infinite1(i,j,i_elmt) /= 0.0_CUSTOM_REAL) then
          !call MatSetValues(Amat1,1,ggdof_elmt(i),1,ggdof_elmt(j), &
          !                  storekmat_infinite1(i,j,i_elmt),ADD_VALUES,ierr); CHECK_PETSC_ERROR(ierr)
          ! petsc types
          v = storekmat_infinite1(i,j,i_elmt)
          call MatSetValues(Amat1,1,ggdof_elmt(i),1,ggdof_elmt(j),v,ADD_VALUES,ierr); CHECK_PETSC_ERROR(ierr)
        endif
      enddo
    enddo
  enddo
  call synchronize_all()
  deallocate(storekmat_infinite1)

  ! user output
  if (myrank == 0) then
    tCPU = wtime() - tstart
    write(IMAIN,*) '    elapsed time: ',sngl(tCPU),'(s)'
    write(IMAIN,*) '    assembling matrix'
    call flush_IMAIN()
  endif

  call MatAssemblyBegin(Amat1,MAT_FINAL_ASSEMBLY,ierr); CHECK_PETSC_ERROR(ierr)
  call MatAssemblyEnd(Amat1,MAT_FINAL_ASSEMBLY,ierr); CHECK_PETSC_ERROR(ierr)
  ! symmetric
  call MatSetOption(Amat1,MAT_SYMMETRIC,PETSC_TRUE,ierr); CHECK_PETSC_ERROR(ierr)

  ! debugging output
  if (DEBUG_FILE_OUTPUT) then
    allocate(cols(nzeros_max))
    cols(:) = 0
    write(char_myrank,'(i4)') myrank
    outf_name='tmp_nonzeros'//trim(adjustl(char_myrank))
    open(1,file=outf_name,action='write',status='replace')
    call MatGetOwnershipRange(Amat1,istart,iend,ierr); CHECK_PETSC_ERROR(ierr)
    do i = istart,iend-1
      cols(:) = -1
      ! gets row i
      call MatGetRow(Amat1,i,ncols,cols,PETSC_NULL_SCALAR,ierr); CHECK_PETSC_ERROR(ierr)
      ndiag = count(cols >= ig0 .and. cols <= ig1)
      noffdiag = ncols-ndiag
      write(1,*) ndiag,noffdiag,ncols
      ! free temporary space of MatGetRow()
      call MatRestoreRow(Amat1,i,ncols,PETSC_NULL_INTEGER,PETSC_NULL_SCALAR,ierr); CHECK_PETSC_ERROR(ierr)
    enddo
    close(1)
    deallocate(cols)
  endif

  ! synchronize all processes
  call synchronize_all()

  ! matrix info
  call MatGetInfo(Amat1, MAT_GLOBAL_MAX, info, ierr); CHECK_PETSC_ERROR(ierr)

  mallocsval = info(MAT_INFO_MALLOCS)               ! number of mallocs during MatSetValues()
  !memval = info(MAT_INFO_MEMORY)                   ! memory allocated - not provided
  !nonzeros_allocated = info(MAT_INFO_NZ_ALLOCATED) ! nonzero entries allocated

  ! memory usage (returns bytes)
  call PetscMemoryGetCurrentUsage(bytes, ierr); CHECK_PETSC_ERROR(ierr)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    number of mallocs during setting values: ',int(mallocsval)
    write(IMAIN,*) '    current PETSc memory usage  = ',sngl(bytes / 1024.d0 / 1024.d0),'MB'
    write(IMAIN,*) '                                = ',sngl(bytes / 1024.d0 / 1024.d0 / 1024.d0),'GB'
    write(IMAIN,*)
    tCPU = wtime() - tstart
    write(IMAIN,*) '    Elapsed time for PETSc solver matrix setup: ',sngl(tCPU),'(s)'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

#else
  ! no PETSc compilation support
  ! compilation without PETSc support
  if (myrank == 0) then
    print *, "Error: PETSc solver enabled without PETSc Support."
    print *, "       To enable PETSc support, reconfigure with --with-petsc flag."
  endif
  ! safety stop
  call exit_MPI(myrank,"Error PETSc solver: petsc_set_matrix1 called without compilation support")

#endif

  end subroutine petsc_set_matrix1

!
!===============================================================================
!


  subroutine petsc_set_vector1(rload1)

#ifdef USE_PETSC
  use constants, only: IFLAG_IN_FICTITIOUS_CUBE,NGLLCUBE_INF

  use specfem_par_full_gravity, only: l2gdof1
#endif

  implicit none

  !PetscScalar,intent(in) :: rload1(0:)
  real(kind=CUSTOM_REAL), intent(in) :: rload1(0:)

#ifdef USE_PETSC
  PetscScalar :: zero

  ! types required by VecSetValues:
  !   VecSetValues(Vec x, PetscInt ni, const PetscInt ix[], const PetscScalar y[], InsertMode iora)
  ! the scalar array must be a PetscScalar type, this might differ from a real(kind=CUSTOM_REAL) type
  PetscScalar :: y(0:size(rload1)-1)

  y(0:) = rload1(0:)

  zero = 0.0
  call VecSet(bvec1,zero,ierr)
  !call VecSetValues(bvec1,neq1,l2gdof1(1:),rload1(1:),ADD_VALUES,ierr)
  call VecSetValues(bvec1,neq1,l2gdof1(1:),y(1:),ADD_VALUES,ierr);

  ! assemble vector
  call VecAssemblyBegin(bvec1,ierr)
  call VecAssemblyEnd(bvec1,ierr)

#else
  ! no PETSc compilation support
  integer :: idummy

  ! compilation without PETSc support
  if (myrank == 0) then
    print *, "Error: PETSc solver enabled without PETSc Support."
    print *, "       To enable PETSc support, reconfigure with --with-petsc flag."
  endif
  ! safety stop
  call exit_MPI(myrank,"Error PETSc solver: petsc_set_vector1 called without compilation support")

  ! to avoid compiler warning
  idummy = size(rload1)
#endif

  end subroutine petsc_set_vector1

!
!===============================================================================
!

! not used so far...

!  subroutine petsc_set_backward_vector1(b_rload1)
!
!#ifdef USE_PETSC
!  use specfem_par_full_gravity, only: l2gdof1, neq1
!#endif
!
!  implicit none
!  !PetscScalar,intent(in) :: b_rload1(0:)
!  real(kind=CUSTOM_REAL), intent(in) :: b_rload1(0:)
!
!#ifdef USE_PETSC
!  PetscScalar :: zero
!  ! types required by VecSetValues:
!  !   VecSetValues(Vec x, PetscInt ni, const PetscInt ix[], const PetscScalar y[], InsertMode iora)
!  ! the scalar array must be a PetscScalar type, this might differ from a real(kind=CUSTOM_REAL) type
!  PetscScalar :: y(0:size(b_rload1)-1)
!
!  y(0:) = b_rload1(0:)
!
!  zero = 0.0
!  call VecSet(b_bvec1,zero,ierr)
!  !call VecSetValues(b_bvec1,neq1,l2gdof1(1:),b_rload1(1:),ADD_VALUES,ierr)
!  call VecSetValues(b_bvec1,neq1,l2gdof1(1:),y(1:),ADD_VALUES,ierr);
!
!  ! assemble vector
!  call VecAssemblyBegin(b_bvec1,ierr)
!  call VecAssemblyEnd(b_bvec1,ierr)
!
!#else
!  ! no PETSc compilation support
!  ! compilation without PETSc support
!  if (myrank == 0) then
!    print *, "Error: PETSc solver enabled without PETSc Support."
!    print *, "       To enable PETSc support, reconfigure with --with-petsc flag."
!  endif
!  ! safety stop
!  call exit_MPI(myrank,"Error PETSc solver: petsc_set_backward_vector1 called without compilation support")
!
!#endif
!
!  end subroutine petsc_set_backward_vector1

!
!===============================================================================
!

  subroutine petsc_solve1(sdata1,niter)

  implicit none
  !PetscScalar :: sdata1(:)
  real(kind=CUSTOM_REAL) :: sdata1(:)
  integer, optional :: niter

#ifdef USE_PETSC
  ! local parameters
  PetscInt :: iter,ireason
  ! petsc type array
  PetscScalar :: y(size(sdata1))

  call KSPSolve(ksp1,bvec1,xvec1,ierr); CHECK_PETSC_ERROR(ierr)
  call KSPGetConvergedReason(ksp1,ireason,ierr); CHECK_PETSC_ERROR(ierr)
  call KSPGetIterationNumber(ksp1,iter,ierr); CHECK_PETSC_ERROR(ierr)

  !debug
  !if (myrank == 0) print *,'debug: petsc_solve1: converged reason: ',ireason

  !call scatter_globalvec1_backward(xvec1, sdata1)

  ! explict conversion to PetscScalar array
  y(:) = sdata1(:)

  call scatter_globalvec1(xvec1, y)

  ! return values
  sdata1(:) = real(y(:),kind=CUSTOM_REAL)

  ! returns number of iterations
  if (present(niter)) niter = iter

#else
  ! no PETSc compilation support
  integer :: idummy

  ! compilation without PETSc support
  if (myrank == 0) then
    print *, "Error: PETSc solver enabled without PETSc Support."
    print *, "       To enable PETSc support, reconfigure with --with-petsc flag."
  endif
  ! safety stop
  call exit_MPI(myrank,"Error PETSc solver: petsc_solve1 called without compilation support")

  ! to avoid compiler warning
  idummy = size(sdata1)
  if (present(niter)) niter = 0
#endif

  end subroutine petsc_solve1

!
!===============================================================================
!

! not used so far...

!  subroutine petsc_backward_solve1(b_sdata1)
!
!  implicit none
!  !PetscScalar :: b_sdata1(:)
!  real(kind=CUSTOM_REAL) :: b_sdata1(:)
!
!#ifdef USE_PETSC
!  ! local parameters
!  PetscInt :: b_iter,b_ireason
!  ! petsc type array
!  PetscScalar :: y(size(b_sdata1))
!
!  call KSPSolve(b_ksp1,b_bvec1,b_xvec1,ierr); CHECK_PETSC_ERROR(ierr)
!  call KSPGetConvergedReason(b_ksp1,b_ireason,ierr); CHECK_PETSC_ERROR(ierr)
!  call KSPGetIterationNumber(b_ksp1,b_iter,ierr); CHECK_PETSC_ERROR(ierr)
!
!  !debug
!  !if (myrank == 0) print *,'debug: petsc_backward_solve1: converged reason: ',b_ireason
!
!  !call scatter_globalvec1_backward(b_xvec1, b_sdata1)
!
!  ! explict conversion to PetscScalar array
!  y(:) = b_sdata1(:)
!
!  call scatter_globalvec1_backward(b_xvec1, y)
!
!  ! return values
!  b_sdata1(:) = y(:)
!
!#else
!  ! no PETSc compilation support
!  ! compilation without PETSc support
!  if (myrank == 0) then
!    print *, "Error: PETSc solver enabled without PETSc Support."
!    print *, "       To enable PETSc support, reconfigure with --with-petsc flag."
!  endif
!  ! safety stop
!  call exit_MPI(myrank,"Error PETSc solver: petsc_backward_solve1 called without compilation support")
!
!#endif
!
!  end subroutine petsc_backward_solve1

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
!  call VecScatterBegin(l_vscat1,global_vec,l_vec1,INSERT_VALUES,SCATTER_FORWARD,ierr); CHECK_PETSC_ERROR(ierr)
!  call VecScatterEnd(l_vscat1,global_vec,l_vec1,INSERT_VALUES,SCATTER_FORWARD,ierr); CHECK_PETSC_ERROR(ierr)
!  call VecGetSize(l_vec1,n,ierr)
!  call VecGetArrayF90(l_vec1,array_data,ierr); CHECK_PETSC_ERROR(ierr)
!  larray(1:n)=array_data(1:n)
!  call VecRestoreArrayF90(l_vec1,array_data,ierr); CHECK_PETSC_ERROR(ierr)
!  return
!  end subroutine scatter_globalvec1

!
!===============================================================================
!

#ifdef USE_PETSC

  subroutine scatter_globalvec1(global_vec,larray)

  implicit none

  type(tVec) :: global_vec
  PetscScalar :: larray(:)

  PetscInt :: n
  PetscScalar,pointer :: array_data(:)

  call VecScatterBegin(vscat1,global_vec,local_vec1,INSERT_VALUES,SCATTER_FORWARD,ierr); CHECK_PETSC_ERROR(ierr)
  call VecScatterEnd(vscat1,global_vec,local_vec1,INSERT_VALUES,SCATTER_FORWARD,ierr); CHECK_PETSC_ERROR(ierr)
  call VecGetSize(local_vec1,n,ierr)

  call VecGetArrayF90(local_vec1,array_data,ierr); CHECK_PETSC_ERROR(ierr)

  larray(1:n) = array_data(1:n)
  call VecRestoreArrayF90(local_vec1,array_data,ierr); CHECK_PETSC_ERROR(ierr)

  end subroutine scatter_globalvec1

#endif

!
!===============================================================================
!

#ifdef USE_PETSC

! not used so far...

!  subroutine scatter_globalvec1_backward(b_global_vec,b_larray)
!
!  implicit none
!
!  type(tVec) :: b_global_vec
!  PetscScalar :: b_larray(:)
!
!  PetscInt :: b_n
!  PetscScalar,pointer :: b_array_data(:)
!
!  call VecScatterBegin(b_vscat1, b_global_vec, b_local_vec1,INSERT_VALUES,SCATTER_FORWARD,ierr); CHECK_PETSC_ERROR(ierr)
!  call VecScatterEnd(b_vscat1, b_global_vec, b_local_vec1,INSERT_VALUES,SCATTER_FORWARD,ierr); CHECK_PETSC_ERROR(ierr)
!  call VecGetSize(b_local_vec1, b_n,ierr)
!
!  call VecGetArrayF90(b_local_vec1, b_array_data,ierr); CHECK_PETSC_ERROR(ierr)
!
!  b_larray(1:b_n) = b_array_data(1:b_n)
!  call VecRestoreArrayF90(b_local_vec1, b_array_data,ierr); CHECK_PETSC_ERROR(ierr)
!
!  end subroutine scatter_globalvec1_backward

#endif

!
!===============================================================================
!

  subroutine petsc_zero_initialguess1()

#ifdef USE_PETSC
  implicit none
  PetscScalar :: zero

  zero = 0.0
  call VecSet(xvec1,zero,ierr)

  ! assemble vector
  call VecAssemblyBegin(xvec1,ierr)
  call VecAssemblyEnd(xvec1,ierr)

#else
  ! no PETSc compilation support
  ! compilation without PETSc support
  if (myrank == 0) then
    print *, "Error: PETSc solver enabled without PETSc Support."
    print *, "       To enable PETSc support, reconfigure with --with-petsc flag."
  endif
  ! safety stop
  call exit_MPI(myrank,"Error PETSc solver: petsc_zero_initialguess1 called without compilation support")

#endif

  end subroutine petsc_zero_initialguess1

!
!===============================================================================
!

  subroutine petsc_zero_backwards_initialguess1()

#ifdef USE_PETSC
  implicit none
  PetscScalar :: zero

  zero = 0.0
  call VecSet(b_xvec1,zero,ierr)

  ! assemble vector
  call VecAssemblyBegin(b_xvec1,ierr)
  call VecAssemblyEnd(b_xvec1,ierr)

#else
  ! no PETSc compilation support
  ! compilation without PETSc support
  if (myrank == 0) then
    print *, "Error: PETSc solver enabled without PETSc Support."
    print *, "       To enable PETSc support, reconfigure with --with-petsc flag."
  endif
  ! safety stop
  call exit_MPI(myrank,"Error PETSc solver: petsc_zero_backwards_initialguess1 called without compilation support")

#endif

  end subroutine petsc_zero_backwards_initialguess1

!
!===============================================================================
!

!  subroutine petsc_set_initialguess1(loc_pgrav1)
!  implicit none
!  PetscScalar,intent(in) :: loc_pgrav1(0:)
!  PetscScalar :: zero
!
!  zero = 0.0
!  call VecSet(bvec1,zero,ierr)
!  call VecSetValues(bvec1,neq1,l2gdof1(1:),loc_pgrav1(1:),INSERT_VALUES,ierr)
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
!  implicit none
!  PetscScalar,intent(in) :: loc_pgrav1(0:)
!  PetscScalar :: zero
!
!  zero = 0.0
!  call VecSet(b_bvec1,zero,ierr)
!  call VecSetValues(b_bvec1, neq1, l2gdof1(1:), loc_pgrav1(1:), INSERT_VALUES, ierr)
!
!  ! assemble vector
!  call VecAssemblyBegin(b_bvec1,ierr)
!  call VecAssemblyEnd(b_bvec1,ierr)
!
!  end subroutine petsc_set_backward_initialguess1

!
!===============================================================================
!

  subroutine petsc_finalize1()

#ifdef USE_PETSC
  use specfem_par, only: SIMULATION_TYPE,USE_POISSON_SOLVER_5GLL

  implicit none

  ! Free work space.  All PETSc objects should be destroyed when they
  ! are no longer needed.

  call VecDestroy(xvec1,ierr); CHECK_PETSC_ERROR(ierr)
  call VecDestroy(uvec1,ierr); CHECK_PETSC_ERROR(ierr)
  call VecDestroy(bvec1,ierr); CHECK_PETSC_ERROR(ierr)
  call MatDestroy(Amat1,ierr); CHECK_PETSC_ERROR(ierr)
  call KSPDestroy(ksp1,ierr); CHECK_PETSC_ERROR(ierr)
  call VecScatterDestroy(vscat1,ierr); CHECK_PETSC_ERROR(ierr)

  if (SIMULATION_TYPE == 3) then
    call VecDestroy(b_xvec1,ierr); CHECK_PETSC_ERROR(ierr)
    call VecDestroy(b_bvec1,ierr); CHECK_PETSC_ERROR(ierr)
    call KSPDestroy(b_ksp1,ierr); CHECK_PETSC_ERROR(ierr)
    call VecScatterDestroy(b_vscat1,ierr); CHECK_PETSC_ERROR(ierr)
  endif

  ! final petsc - otherwise called in petsc_finalize() routine
  if (.not. USE_POISSON_SOLVER_5GLL) then
    call PetscFinalize(ierr); CHECK_PETSC_ERROR(ierr)
  endif

#else
  ! no PETSc compilation support
  ! compilation without PETSc support
  if (myrank == 0) then
    print *, "Error: PETSc solver enabled without PETSc Support."
    print *, "       To enable PETSc support, reconfigure with --with-petsc flag."
  endif
  ! safety stop
  call exit_MPI(myrank,"Error PETSc solver: petsc_finalize1 called without compilation support")

#endif

  end subroutine petsc_finalize1


!===============================================================================
! Level-2 solver
!===============================================================================

  subroutine petsc_initialize()

#ifdef USE_PETSC
  implicit none
  PetscInt :: istart,iend
  PetscInt, allocatable :: nzeros(:)
  ! memory info
  PetscLogDouble :: bytes

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    initializing PETSc Level-2 solver'
    call flush_IMAIN()
  endif

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Beginning of program
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr); CHECK_PETSC_ERROR(ierr)

  !call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-n',ngdof,flg,ierr)
  !call synchronize_all()

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Compute the matrix and right-hand-side vector that define
  ! the linear system, Ax = b.
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ! Create matrix. When using MatCreate(), the matrix format can
  ! be specified at runtime.

  ! count number of nonzeros per row
  allocate(nzeros(neq))
  nzeros(:) = 0
  nzeros(krow_sparse(:)) = nzeros(krow_sparse(:)) + 1

  nzeros_max = maxvec(nzeros)
  nzeros_min = minvec(nzeros)

  !debug
  if (myrank == 0) print *,'PETSc Level-2 solver: ngdof:',ngdof,' nzeros_max:',nzeros_max,' nzeros_min:',nzeros_min

  call MatCreate(PETSC_COMM_WORLD,Amat,ierr); CHECK_PETSC_ERROR(ierr)
  call MatSetType(Amat,MATMPIAIJ,ierr); CHECK_PETSC_ERROR(ierr)
  call MatSetSizes(Amat,PETSC_DECIDE,PETSC_DECIDE,ngdof,ngdof,ierr); CHECK_PETSC_ERROR(ierr)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    number of global degrees of freedom   : ',ngdof
    write(IMAIN,*) '    number of non-zero row entries min/max: ',nzeros_min,'/',nzeros_max
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! preallocation
  call MatMPIAIJSetPreallocation(Amat,nzeros_max,nzeros,nzeros_max,20*nzeros,ierr); CHECK_PETSC_ERROR(ierr)
  call MatSetFromOptions(Amat,ierr); CHECK_PETSC_ERROR(ierr)
  call MatGetOwnershipRange(Amat,istart,iend,ierr); CHECK_PETSC_ERROR(ierr)

  !debug
  if (myrank == 0) print *,'PETSc Level-2 solver: actual global index range:',minval(kgrow_sparse),maxval(kgrow_sparse)

  !debug
  print *,'PETSc Level-2 solver: global index:',myrank,istart,iend,iend-istart
  call synchronize_all()

  ! free temporary array
  deallocate(nzeros)
  deallocate(krow_sparse,kcol_sparse)
  deallocate(kgrow_sparse,kgcol_sparse)

  !debug
  if (myrank == 0) print *,'PETSc Level-2 solver: matrix'
  call synchronize_all()

  ! Create vectors.  Note that we form 1 vector from scratch and
  ! then duplicate as needed.

  !call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,ngdof,xvec,ierr)
  call VecCreate(PETSC_COMM_WORLD,xvec,ierr); CHECK_PETSC_ERROR(ierr)
  call VecSetSizes(xvec,PETSC_DECIDE,ngdof,ierr); CHECK_PETSC_ERROR(ierr)
  call VecSetFromOptions(xvec,ierr); CHECK_PETSC_ERROR(ierr)
  call VecDuplicate(xvec,bvec,ierr); CHECK_PETSC_ERROR(ierr)
  call VecDuplicate(xvec,uvec,ierr); CHECK_PETSC_ERROR(ierr)

  !debug
  if (myrank == 0) print *,'PETSc Level-2 solver: vector'

  ! memory usage (returns bytes)
  call PetscMemoryGetCurrentUsage(bytes, ierr); CHECK_PETSC_ERROR(ierr)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    current PETSc memory usage  = ',sngl(bytes / 1024.d0 / 1024.d0),'MB'
    write(IMAIN,*) '                                = ',sngl(bytes / 1024.d0 / 1024.d0 / 1024.d0),'GB'
    call flush_IMAIN()
  endif

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Create the linear solver and set various options
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ! Create linear solver context

  call KSPCreate(PETSC_COMM_WORLD,ksp,ierr); CHECK_PETSC_ERROR(ierr)

  ! Set operators. Here the matrix that defines the linear system
  ! also serves as the preconditioning matrix.
  !call KSPSetOperators(ksp,Amat,Amat,SAME_PRECONDITIONER,ierr); CHECK_PETSC_ERROR(ierr) ! version < 3.5
  call KSPSetOperators(ksp,Amat,Amat,ierr); CHECK_PETSC_ERROR(ierr) ! version >= 3.5

  call KSPSetType(ksp,KSPCG,ierr); CHECK_PETSC_ERROR(ierr)

  !debug
  if (myrank == 0) print *,'PETSc Level-2 solver: ksp0'

  call KSPGetPC(ksp,pc,ierr); CHECK_PETSC_ERROR(ierr)
  call PCSetType(pc,PCHYPRE,ierr); CHECK_PETSC_ERROR(ierr)

  !debug
  if (myrank == 0) print *,'PETSc Level-2 solver: ksp1'

  call KSPSetTolerances(ksp,KSP_RTOL,KSP_ATOL,KSP_DTOL,KSP_MAXITER,ierr); CHECK_PETSC_ERROR(ierr)

  !debug
  if (myrank == 0) print *,'PETSc Level-2 solver: ksp2'

  !  Set runtime options, e.g.,
  !    -ksp_type < type> -pc_type < type> -ksp_monitor -ksp_KSP_RTOL < KSP_RTOL>
  !  These options will override those specified above as long as
  !  KSPSetFromOptions() is called _after_ any other customization
  !  routines.
  call KSPSetFromOptions(ksp,ierr); CHECK_PETSC_ERROR(ierr)

  !debug
  if (myrank == 0) print *,'PETSc Level-2 solver: initialization done.'

  ! synchronize all processes
  call synchronize_all()

#else
  ! no PETSc compilation support
  ! compilation without PETSc support
  if (myrank == 0) then
    print *, "Error: PETSc solver enabled without PETSc Support."
    print *, "       To enable PETSc support, reconfigure with --with-petsc flag."
  endif
  ! safety stop
  call exit_MPI(myrank,"Error PETSc solver: petsc_initialize called without compilation support")

#endif

  end subroutine petsc_initialize

!
!===============================================================================
!

  subroutine petsc_set_matrix()

#ifdef USE_PETSC
  use constants, only: CUSTOM_REAL,IFLAG_IN_FICTITIOUS_CUBE,NEDOF

  use specfem_par, only: NSPEC_INNER_CORE, &
    NSPEC_OUTER_CORE,NSPEC_CRUST_MANTLE,NSPEC_TRINFINITE,NSPEC_INFINITE

  use specfem_par_innercore, only: idoubling_inner_core

  use specfem_par_full_gravity, only: &
    ggdof_oc,storekmat_outer_core,inode_elmt_oc, &
    ggdof_ic,storekmat_inner_core,inode_elmt_ic, &
    ggdof_cm,storekmat_crust_mantle,inode_elmt_cm, &
    ggdof_trinf,storekmat_trinfinite,inode_elmt_trinf, &
    ggdof_inf,storekmat_infinite,inode_elmt_inf

  implicit none
  integer :: i,i_elmt,j,ncount
  integer :: ggdof_elmt(NEDOF),idof(NEDOF),igdof(NEDOF)

  ! types required by MatSetValues:
  !   MatSetValues(Mat mat, PetscInt m, const PetscInt idxm[], PetscInt n, \
  !                const PetscInt idxn[], const PetscScalar v[], InsertMode addv)
  PetscScalar :: v
  PetscScalar,dimension(:,:),allocatable :: varr

  ! matrix info
  double precision :: info(MAT_INFO_SIZE)
  double precision :: mallocsval
  PetscLogDouble :: bytes

  ! timing
  double precision :: tstart,tCPU
  double precision, external :: wtime

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    setting up Level-2 solver matrix'
    call flush_IMAIN()
  endif

  ! timing
  tstart = wtime()

  ! Set and assemble matrix.
  !  - Note that MatSetValues() uses 0-based row and column numbers
  !  in Fortran as well as in C (as set here in the array "col").
    ! stage 0: store all elements

  call MatZeroEntries(Amat,ierr); CHECK_PETSC_ERROR(ierr)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    setup inner core values'
    call flush_IMAIN()
  endif

  ! inner core
  do i_elmt = 1,NSPEC_INNER_CORE
    if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE)cycle
    ggdof_elmt(:) = reshape(ggdof_ic(:,inode_elmt_ic(:,i_elmt)),(/NEDOF/))
    ggdof_elmt(:) = ggdof_elmt(:)-1 ! petsc index starts from 0
    !ncount = 0; idof = -1; igdof = -1
    do i = 1,NEDOF
      do j = 1,NEDOF
        if (ggdof_elmt(i) >= 0 .and. ggdof_elmt(j) >= 0 .and. &
            storekmat_inner_core(i,j,i_elmt) /= 0.0_CUSTOM_REAL) then
          ! this leads to an error of no specific subroutine found
          !call MatSetValues(Amat,1,ggdof_elmt(i),1,ggdof_elmt(j), &
          !                  storekmat_inner_core(i,j,i_elmt),ADD_VALUES,ierr); CHECK_PETSC_ERROR(ierr)
          ! petsc types
          v = storekmat_inner_core(i,j,i_elmt)
          call MatSetValues(Amat,1,ggdof_elmt(i),1,ggdof_elmt(j),v,ADD_VALUES,ierr); CHECK_PETSC_ERROR(ierr)
        endif
      enddo
    enddo
  enddo
  deallocate(storekmat_inner_core)

  ! user output
  if (myrank == 0) then
    tCPU = wtime() - tstart
    write(IMAIN,*) '    elapsed time: ',sngl(tCPU),'(s)'
    write(IMAIN,*) '    setup outer core values'
    call flush_IMAIN()
  endif

  ! outer core
  do i_elmt = 1,NSPEC_OUTER_CORE
    ggdof_elmt(:) = reshape(ggdof_oc(:,inode_elmt_oc(:,i_elmt)),(/NEDOF/))
    ggdof_elmt(:) = ggdof_elmt(:)-1 ! petsc index starts from 0
    !ncount = 0; idof = -1; igdof = -1
    do i = 1,NEDOF
      do j = 1,NEDOF
        if (ggdof_elmt(i) >= 0 .and. ggdof_elmt(j) >= 0 .and. &
            storekmat_outer_core(i,j,i_elmt) /= 0.0_CUSTOM_REAL) then
          ! petsc types
          v = storekmat_outer_core(i,j,i_elmt)
          call MatSetValues(Amat,1,ggdof_elmt(i),1,ggdof_elmt(j),v,ADD_VALUES,ierr); CHECK_PETSC_ERROR(ierr)
        endif
      enddo
    enddo
  enddo
  deallocate(storekmat_outer_core)

  ! user output
  if (myrank == 0) then
    tCPU = wtime() - tstart
    write(IMAIN,*) '    elapsed time: ',sngl(tCPU),'(s)'
    write(IMAIN,*) '    setup crust/mantle values'
    call flush_IMAIN()
  endif

  ! crust mantle
  do i_elmt = 1,NSPEC_CRUST_MANTLE
    ggdof_elmt(:) = reshape(ggdof_cm(:,inode_elmt_cm(:,i_elmt)),(/NEDOF/))
    ggdof_elmt(:) = ggdof_elmt(:)-1 ! petsc index starts from 0
    ncount = 0; idof(:) = -1; igdof(:) = -1
    do i = 1,NEDOF
      if (ggdof_elmt(i) >= 0) then
        ncount = ncount+1
        idof(ncount) = i
        igdof(ncount) = ggdof_elmt(i)
      endif
    enddo
    !call MatSetValues(Amat,ncount,igdof(1:ncount),ncount,igdof(1:ncount), &
    !                  storekmat_crust_mantle(idof(1:ncount),idof(1:ncount),i_elmt),ADD_VALUES,ierr); CHECK_PETSC_ERROR(ierr)
    ! petsc types
    allocate(varr(ncount,ncount))
    varr(:,:) = storekmat_crust_mantle(idof(1:ncount),idof(1:ncount),i_elmt)
    call MatSetValues(Amat,ncount,igdof(1:ncount),ncount,igdof(1:ncount),varr,ADD_VALUES,ierr); CHECK_PETSC_ERROR(ierr)
    deallocate(varr)
  enddo
  deallocate(storekmat_crust_mantle)

  ! user output
  if (myrank == 0) then
    tCPU = wtime() - tstart
    write(IMAIN,*) '    elapsed time: ',sngl(tCPU),'(s)'
    call flush_IMAIN()
  endif

  ! trinfinite
  if (NSPEC_TRINFINITE > 0) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '    setup trinfinite values'
      call flush_IMAIN()
    endif

    do i_elmt = 1,NSPEC_TRINFINITE
      ggdof_elmt(:) = reshape(ggdof_trinf(:,inode_elmt_trinf(:,i_elmt)),(/NEDOF/))
      ggdof_elmt(:) = ggdof_elmt-1 ! petsc index starts from 0
      ncount = 0; idof(:) = -1; igdof(:) = -1
      do i = 1,NEDOF
        if (ggdof_elmt(i) >= 0) then
          ncount = ncount+1
          idof(ncount) = i
          igdof(ncount) = ggdof_elmt(i)
        endif
      enddo
      !call MatSetValues(Amat,ncount,igdof(1:ncount),ncount,igdof(1:ncount), &
      !                  storekmat_trinfinite(idof(1:ncount),idof(1:ncount),i_elmt),ADD_VALUES,ierr); CHECK_PETSC_ERROR(ierr)
      ! petsc types
      allocate(varr(ncount,ncount))
      varr(:,:) = storekmat_trinfinite(idof(1:ncount),idof(1:ncount),i_elmt)
      call MatSetValues(Amat,ncount,igdof(1:ncount),ncount,igdof(1:ncount),varr,ADD_VALUES,ierr); CHECK_PETSC_ERROR(ierr)
      deallocate(varr)
    enddo
    ! user output
    if (myrank == 0) then
      tCPU = wtime() - tstart
      write(IMAIN,*) '    elapsed time: ',sngl(tCPU),'(s)'
      call flush_IMAIN()
    endif
  endif
  deallocate(storekmat_trinfinite)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    setup infinite values'
    call flush_IMAIN()
  endif

  ! infinite
  do i_elmt = 1,NSPEC_INFINITE
    ggdof_elmt(:) = reshape(ggdof_inf(:,inode_elmt_inf(:,i_elmt)),(/NEDOF/))
    ggdof_elmt(:) = ggdof_elmt(:)-1 ! petsc index starts from 0
    ncount = 0; idof(:) = -1; igdof(:) = -1
    do i = 1,NEDOF
      if (ggdof_elmt(i) >= 0) then
        ncount = ncount+1
        idof(ncount) = i
        igdof(ncount) = ggdof_elmt(i)
      endif
    enddo
    !call MatSetValues(Amat,ncount,igdof(1:ncount),ncount,igdof(1:ncount), &
    !                  storekmat_infinite(idof(1:ncount),idof(1:ncount),i_elmt),ADD_VALUES,ierr); CHECK_PETSC_ERROR(ierr)
    ! petsc types
    allocate(varr(ncount,ncount))
    varr(:,:) = storekmat_infinite(idof(1:ncount),idof(1:ncount),i_elmt)
    call MatSetValues(Amat,ncount,igdof(1:ncount),ncount,igdof(1:ncount),varr,ADD_VALUES,ierr); CHECK_PETSC_ERROR(ierr)
    deallocate(varr)
  enddo
  deallocate(storekmat_infinite)

  ! user output
  if (myrank == 0) then
    tCPU = wtime() - tstart
    write(IMAIN,*) '    elapsed time: ',sngl(tCPU),'(s)'
    write(IMAIN,*) '    assembling matrix'
    call flush_IMAIN()
  endif

  call MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY,ierr); CHECK_PETSC_ERROR(ierr)
  call MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY,ierr); CHECK_PETSC_ERROR(ierr)

  ! matrix info
  call MatGetInfo(Amat, MAT_GLOBAL_MAX, info, ierr); CHECK_PETSC_ERROR(ierr)
  mallocsval = info(MAT_INFO_MALLOCS) ! number of mallocs during MatSetValues()

  ! memory usage
  call PetscMemoryGetCurrentUsage(bytes, ierr); CHECK_PETSC_ERROR(ierr)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    number of mallocs during setting values: ',int(mallocsval)
    write(IMAIN,*) '    current PETSc memory usage  = ',sngl(bytes / 1024.d0 / 1024.d0),'MB'
    write(IMAIN,*) '                                = ',sngl(bytes / 1024.d0 / 1024.d0 / 1024.d0),'GB'
    write(IMAIN,*)
    tCPU = wtime() - tstart
    write(IMAIN,*) '    Elapsed time for PETSc solver matrix setup: ',sngl(tCPU),'(s)'
    write(IMAIN,*)
    call flush_IMAIN()

    call flush_IMAIN()
  endif

#else
  ! no PETSc compilation support
  ! compilation without PETSc support
  if (myrank == 0) then
    print *, "Error: PETSc solver enabled without PETSc Support."
    print *, "       To enable PETSc support, reconfigure with --with-petsc flag."
  endif
  ! safety stop
  call exit_MPI(myrank,"Error PETSc solver: petsc_set_matrix called without compilation support")

#endif

  end subroutine petsc_set_matrix

!
!===============================================================================
!

  subroutine petsc_set_vector()

#ifdef USE_PETSC
  use specfem_par_full_gravity, only: l2gdof,gravload

  implicit none
  PetscScalar :: zero

  ! petsc explicit types
  PetscInt :: ni
  PetscInt :: ix(size(l2gdof(1:)))
  PetscScalar :: y(size(gravload(1:)))

  ! Set exact solution; then compute right-hand-side vector.
  !none=-1.0
  !one=1.0
  zero = 0.0
  call VecSet(bvec,zero,ierr); CHECK_PETSC_ERROR(ierr)

  ! types required by VecSetValues:
  !   VecSetValues(Vec x, PetscInt ni, const PetscInt ix[], const PetscScalar y[], InsertMode iora)
  ! this leads to an error of no specific subroutine found
  !call VecSetValues(bvec,neq,l2gdof(1:),gravload(1:),ADD_VALUES,ierr); CHECK_PETSC_ERROR(ierr)

  ! sets up petsc types
  ni = neq
  ix(:) = l2gdof(1:)
  y(:) = gravload(1:)

  call VecSetValues(bvec,neq,ix,y,ADD_VALUES,ierr); CHECK_PETSC_ERROR(ierr)

  ! assemble vector
  call VecAssemblyBegin(bvec,ierr); CHECK_PETSC_ERROR(ierr)
  call VecAssemblyEnd(bvec,ierr); CHECK_PETSC_ERROR(ierr)

#else
  ! no PETSc compilation support
  ! compilation without PETSc support
  if (myrank == 0) then
    print *, "Error: PETSc solver enabled without PETSc Support."
    print *, "       To enable PETSc support, reconfigure with --with-petsc flag."
  endif
  ! safety stop
  call exit_MPI(myrank,"Error PETSc solver: petsc_set_vector called without compilation support")

#endif

  end subroutine petsc_set_vector

!
!===============================================================================
!

  subroutine petsc_solve(sdata,niter)

  implicit none
  !PetscScalar :: sdata(:)
  real(kind=CUSTOM_REAL) :: sdata(:)
  integer, optional :: niter

#ifdef USE_PETSC
  ! local parameters
  PetscInt :: iter,ireason
  ! petsc type array
  PetscScalar :: y(size(sdata))

  call KSPSolve(ksp,bvec,xvec,ierr); CHECK_PETSC_ERROR(ierr)

  ! View solver info; we could instead use the option -ksp_view
  call KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr); CHECK_PETSC_ERROR(ierr)

  !-------------------------------------------------------------------------------
  ! Check solution and clean up
  !-------------------------------------------------------------------------------

  call KSPGetConvergedReason(ksp,ireason,ierr); CHECK_PETSC_ERROR(ierr)
  call KSPGetIterationNumber(ksp,iter,ierr); CHECK_PETSC_ERROR(ierr)

  !debug
  if (myrank == 0) then
    print *,'PETSc solve: Converged reason: ',ireason
    print *,'PETSc solve: Iterations      : ',iter
  endif

  ! Check the error
  !call VecAXPY(xvec,none,uvec,ierr); CHECK_PETSC_ERROR(ierr)
  !call VecNorm(xvec,NORM_2,norm,ierr); CHECK_PETSC_ERROR(ierr)
  !
  !if (norm > 1.e-12) then
  !  write(*,'(a,e11.4,a,i5)')'Norm of error:',norm,', Iterations:',its
  !else
  !  write(*,'(a,i5,a)')'Norm of error < 1.e-12, Iterations:',its
  !endif

  ! explict conversion to PetscScalar array
  y(:) = sdata(:)
  ! no scatter done
  !call scatter_globalvec1(xvec1, y)
  ! return values
  !sdata(:) = real(y(:),kind=CUSTOM_REAL)

  ! returns number of iterations
  if (present(niter)) niter = iter

#else
  ! no PETSc compilation support
  integer :: idummy

  ! compilation without PETSc support
  if (myrank == 0) then
    print *, "Error: PETSc solver enabled without PETSc Support."
    print *, "       To enable PETSc support, reconfigure with --with-petsc flag."
  endif
  ! safety stop
  call exit_MPI(myrank,"Error PETSc solver: petsc_solve called without compilation support")

  ! to avoid compiler warning
  idummy = size(sdata)
  if (present(niter)) niter = 0
#endif

  end subroutine petsc_solve

!
!===============================================================================
!

  subroutine petsc_finalize()

#ifdef USE_PETSC
  implicit none

  ! Free work space.  All PETSc objects should be destroyed when they
  ! are no longer needed.

  call VecDestroy(xvec,ierr); CHECK_PETSC_ERROR(ierr)
  call VecDestroy(uvec,ierr); CHECK_PETSC_ERROR(ierr)
  call VecDestroy(bvec,ierr); CHECK_PETSC_ERROR(ierr)
  call MatDestroy(Amat,ierr); CHECK_PETSC_ERROR(ierr)
  call KSPDestroy(ksp,ierr); CHECK_PETSC_ERROR(ierr)

  call PetscFinalize(ierr); CHECK_PETSC_ERROR(ierr)

#else
  ! no PETSc compilation support
  ! compilation without PETSc support
  if (myrank == 0) then
    print *, "Error: PETSc solver enabled without PETSc Support."
    print *, "       To enable PETSc support, reconfigure with --with-petsc flag."
  endif
  ! safety stop
  call exit_MPI(myrank,"Error PETSc solver: petsc_finalize called without compilation support")

#endif

  end subroutine petsc_finalize

!
!===============================================================================
!

#ifdef USE_PETSC

  subroutine check_petsc_err(ierr,line)

  use specfem_par, only: myrank

  implicit none
  PetscErrorCode, intent(in) :: ierr
  integer, intent(in) :: line

  ! the PETSc base type file petscsys.h defines an error checking macro CHKERRA(ierr) as
  ! long free-from version:
  !#define CHKERRA(ierr) if (ierr /= 0) then;call PetscErrorF(ierr,__LINE__,__FILE__);\
  !                      call MPIU_Abort(PETSC_COMM_SELF,ierr);endif
  ! short version:
  !#define CHKERRA(ierr) if (ierr /= 0) then;call PetscErrorF(ierr);\
  !                      call MPIU_Abort(PETSC_COMM_SELF,ierr);endif
  !
  ! given the __LINE__ and __FILE__ info is not very useful when called in this subroutine, we use the short version here.

  if (ierr /= 0) then
    ! petsc error function
#if defined(PETSC_HAVE_FORTRAN_FREE_LINE_LENGTH_NONE)
    call PetscErrorF(ierr,line,"SIEM_solver_petsc.F90")
#else
    call PetscErrorF(ierr)
#endif

    ! user info
    print *
    print *,'Error: PETSc error occurred for rank ',myrank,' in module SIEM_solver_petsc on line ',line
    print *,'       Please check, aborting now...'
    print *

    ! abort
    call MPIU_Abort(PETSC_COMM_SELF,ierr)
  endif

  end subroutine check_petsc_err

#endif

end module siem_solver_petsc

