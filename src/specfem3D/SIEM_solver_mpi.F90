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



module siem_solver_mpi

  use constants, only: myrank,CUSTOM_REAL

  use siem_math_library_mpi, only: maxvec,dot_product_all_proc

  implicit none

  private

  ! maximum number of iteration for conjugate gradient solver
  integer, parameter :: CG_MAXITER = 10000

  ! relative tolerance for conjugate gradient solver
  real(kind=CUSTOM_REAL), parameter :: CG_TOL  = 1.0e-7_CUSTOM_REAL, &
                                       CG_TOL1 = 1.0e-7_CUSTOM_REAL

  ! public functions
  public :: cg_solver
  public :: cg_solver3
  public :: diagpcg_solver
  public :: diagpcg_solver3
  public :: interpolate3to5

contains

!-------------------------------------------------------------------------------

! petsc SP solver preconditioned conjugate-gradient solver

! not used so far...

!  subroutine pcpetsc_cg_solver(myid,neq,u_g,f,dprecon_g)
!
!  implicit none
!  integer,intent(in) :: myid,neq
!  real(kind=CUSTOM_REAL),dimension(0:neq),intent(inout) :: u_g
!  real(kind=CUSTOM_REAL),dimension(0:neq),intent(in) :: f,dprecon_g
!
!  ! local parameters
!  integer :: iter
!  real(kind=CUSTOM_REAL) :: alpha,beta,pkp,rz,maxp,maxu,tol
!  real(kind=CUSTOM_REAL),dimension(0:neq) :: kp,p,p_g,r,r0,z,z_g!,r_g
!
!  real(kind=CUSTOM_REAL),parameter :: zero=0.0_CUSTOM_REAL,zerotol=1.0e-12_CUSTOM_REAL
!
!  ! all global array variables are both MPI and regionally assembled.
!  ! local array variables are regionally assembled.
!  ! for MPI assembly of such array, we have to scatter again to region with
!  ! regionally assembled values.
!
!  ! PCG solver
!
!  ! check if RHS is 0
!  if (maxvec(abs(f)) <= zerotol) then
!    u_g(:) = zero
!    ! all done
!    return
!  endif
!
!  kp = zero
!  if (maxval(abs(u_g)) > zero) then
!    call product_stiffness_vector(neq,u_g,kp)
!  endif
!  ! assemble kp across the regions
!  r = f-kp
!  z = dprecon_g*r
!
!  call scatter_and_assemble(neq,z,z_g)
!
!  p = z
!
!  ! pcg iteration
!  pcg: do iter = 1,CG_MAXITER
!    call scatter_and_assemble(neq,p,p_g)
!
!    call product_stiffness_vector(neq,p_g,kp)
!
!    rz = dot_product_all_proc(r,z_g)
!    pkp = dot_product_all_proc(p_g,kp)
!
!    ! avoids division by zero
!    if (pkp /= zero) then
!      alpha = rz / pkp
!    else
!      alpha = zero
!    endif
!    u_g = u_g + alpha * p_g
!
!    maxp = maxvec(abs(p_g)); maxu = maxvec(abs(u_g))
!
!    ! avoids division by zero
!    if (maxu /= zero) then
!      tol = abs(alpha) * maxp / maxu
!    else
!      tol = abs(alpha) * maxp
!    endif
!
!    if (tol <= CG_TOL) then
!      ! all done
!      return
!    endif
!
!    r0 = r
!    r = r - alpha * kp
!
!    ! solve using single precision petsc solver
!    z = dprecon_g * r
!
!    call scatter_and_assemble(neq,z,z_g)
!
!    !beta = dot_product_all_proc(r,z_g)/rz ! Fletcher–Reeves
!    beta = dot_product_all_proc(r-r0,z_g)/rz !  Polak–Ribière
!
!    p = z + beta * p
!  enddo pcg
!
!  if (myid == 0) print *,'ERROR: PCG solver doesn''t converge! Tolerance:',tol, &
!                         'alpha/maxp/maxu',alpha,maxp,maxu
!  call synchronize_all()
!
!  ! abort
!  call exit_MPI(myrank,'PCG solver does not converge')
!
!  end subroutine pcpetsc_cg_solver

!
!===============================================================================
!

! conjugate-gradient solver

  subroutine cg_solver(myid,neq,u_g,f,niter)

  implicit none
  integer,intent(in) :: myid,neq
  real(kind=CUSTOM_REAL),dimension(0:neq),intent(inout) :: u_g
  real(kind=CUSTOM_REAL),dimension(0:neq),intent(in) :: f
  integer, optional :: niter

  ! local parameters
  integer :: iter
  real(kind=CUSTOM_REAL) :: alpha,beta,pkp,rz,maxp,maxu,tol
  real(kind=CUSTOM_REAL),dimension(0:neq) :: kp,p,p_g
  real(kind=CUSTOM_REAL),dimension(0:neq) :: r,r_g !,r0

  real(kind=CUSTOM_REAL),parameter :: zero=0.0_CUSTOM_REAL,zerotol=1.0e-30_CUSTOM_REAL

  ! all global array variables are both MPI and regionally assembled.
  ! local array variables are regionally assembled.
  ! for MPI assembly of such array, we have to scatter again to region with
  ! regionally assembled values.

  ! PCG solver

  ! check if RHS is 0
  if (maxvec(abs(f)) <= zerotol) then
    u_g(:) = zero
    ! returns number of iterations
    if (present(niter)) niter = 0
    ! all done
    return
  endif

  kp = zero
  if (maxval(abs(u_g)) > zero) then
    call product_stiffness_vector(neq,u_g,kp)
  endif
  ! assemble kp across the regions
  r = f-kp

  call scatter_and_assemble(neq,r,r_g)

  p = r

  ! pcg iteration
  pcg: do iter = 1,CG_MAXITER
    call scatter_and_assemble(neq,p,p_g)

    call product_stiffness_vector(neq,p_g,kp)

    rz = dot_product_all_proc(r,r_g)
    pkp = dot_product_all_proc(p_g,kp)

    ! avoids division by zero
    if (pkp /= zero) then
      alpha = rz / pkp
    else
      alpha = zero
    endif
    u_g = u_g + alpha*p_g

    maxp = maxvec(abs(p_g)); maxu = maxvec(abs(u_g))

    ! avoids division by zero
    if (maxu /= zero) then
      tol = abs(alpha) * maxp / maxu
    else
      tol = abs(alpha) * maxp
    endif

    if (abs(alpha)*maxp/maxu <= CG_TOL) then
      ! returns number of iterations
      if (present(niter)) niter = iter
      ! all done
      return
    endif

    !r0 = r
    r = r - alpha * kp

    call scatter_and_assemble(neq,r,r_g)

    beta = dot_product_all_proc(r,r_g)/rz ! Fletcher–Reeves
    !beta = dot_product_all_proc(r-r0,r_g)/rz !  Polak–Ribière

    p = r + beta * p
  enddo pcg

  if (myid == 0) print *,'ERROR: CG solver doesn''t converge! Tolerance:',tol, &
                         'alpha/maxp/maxu',alpha,maxp,maxu
  call synchronize_all()

  ! abort
  call exit_MPI(myrank,'CG solver does not converge in routine cg_solver()')

  end subroutine cg_solver

!
!============================================
!

! conjugate-gradient solver

  subroutine cg_solver3(myid,neq,u_g,f,niter)

  implicit none
  integer,intent(in) :: myid,neq
  real(kind=CUSTOM_REAL),dimension(0:neq),intent(inout) :: u_g
  real(kind=CUSTOM_REAL),dimension(0:neq),intent(in) :: f
  integer, optional :: niter

  ! local parameters
  integer :: iter
  real(kind=CUSTOM_REAL) :: alpha,beta,pkp,rz,maxf,maxp,maxu,tol
  real(kind=CUSTOM_REAL),dimension(0:neq) :: kp,p,p_g
  real(kind=CUSTOM_REAL),dimension(0:neq) :: r,r_g  !,r0

  real(kind=CUSTOM_REAL),parameter :: zero=0.0_CUSTOM_REAL,zerotol=1.0e-30_CUSTOM_REAL

  ! all global array variables are both MPI and regionally assembled.
  ! local array variables are regionally assembled.
  ! for MPI assembly of such array, we have to scatter again to region with
  ! regionally assembled values.

  ! CG solver
  maxf = maxvec(abs(f))

  !debug
  !if (myid == 0) print *,'CG solver: load: ',maxf

  ! check if RHS is 0
  if (maxf <= zerotol) then
    u_g(:) = zero
    ! returns number of iterations
    if (present(niter)) niter = 0
    ! all done
    return
  endif

  kp = zero
  if (maxval(abs(u_g)) > zero) then
    call product_stiffness_vector3(neq,u_g,kp)
  endif
  ! assemble kp across the regions
  r = f-kp

  call scatter_and_assemble3(neq,r,r_g)

  p = r

  ! pcg iteration
  pcg: do iter = 1,CG_MAXITER
    call scatter_and_assemble3(neq,p,p_g)

    call product_stiffness_vector3(neq,p_g,kp)

    rz = dot_product_all_proc(r,r_g)
    pkp = dot_product_all_proc(p_g,kp)

    !if (abs(pkp) == zero) then
    !  !debug
    !!  if (myid == 0) print *,'CG solver: pkp : ',pkp,' iter: ',iter
    !  ! returns number of iterations
    !  if (present(niter)) niter = iter
    !  ! all done
    !  return
    !endif

    ! avoids division by zero
    if (pkp /= zero) then
      alpha = rz / pkp
    else
      alpha = zero
    endif
    u_g = u_g + alpha * p_g

    maxp = maxvec(abs(p_g)); maxu = maxvec(abs(u_g))

    ! avoids division by zero
    if (maxu /= zero) then
      tol = abs(alpha) * maxp / maxu
    else
      tol = abs(alpha) * maxp
    endif

    ! checks convergence
    if (tol <= CG_TOL1) then
      !debug
      !if (myid == 0) print *,'CG solver: tol : ',tol,' iter: ',iter
      ! returns number of iterations
      if (present(niter)) niter = iter
      ! all done
      return
    endif

    !r0 = r !PR: Polak–Ribière
    r = r - alpha * kp

    call scatter_and_assemble3(neq,r,r_g)

    beta = dot_product_all_proc(r,r_g)/rz ! Fletcher–Reeves
    !beta = dot_product_all_proc(r-r0,r_g)/rz !PR: Polak–Ribière

    p = r + beta * p
  enddo pcg

  if (myid == 0) print *,'ERROR: CG solver doesn''t converge! Tolerance:',tol, &
                         'alpha/maxp/maxu',alpha,maxp,maxu
  call synchronize_all()

  ! abort
  call exit_MPI(myrank,'CG solver does not converge in routine cg_solver3()')

  end subroutine cg_solver3

!
!===============================================================================
!

! diagonally preconditioned conjugate-gradient solver

  subroutine diagpcg_solver(myid,neq,u_g,f,dprecon_g,niter)

  implicit none
  integer,intent(in) :: myid,neq
  real(kind=CUSTOM_REAL),dimension(0:neq),intent(inout) :: u_g
  real(kind=CUSTOM_REAL),dimension(0:neq),intent(in) :: f,dprecon_g
  integer, optional :: niter

  ! local parameters
  integer :: iter
  real(kind=CUSTOM_REAL) :: alpha,beta,pkp,rz,maxp,maxu,tol
  real(kind=CUSTOM_REAL),dimension(0:neq) :: kp,p,p_g,r,r0,z,z_g!,r_g

  real(kind=CUSTOM_REAL),parameter :: zero=0.0_CUSTOM_REAL,zerotol=1.0e-12_CUSTOM_REAL

  ! all global array variables are both MPI and regionally assembled.
  ! local array variables are regionally assembled.
  ! for MPI assembly of such array, we have to scatter again to region with
  ! regionally assembled values.

  ! PCG solver

  ! check if RHS is 0
  if (maxvec(abs(f)) <= zerotol) then
    u_g(:) = zero
    ! returns number of iterations
    if (present(niter)) niter = 0
    ! all done
    return
  endif

  kp = zero
  if (maxval(abs(u_g)) > zero) then
    call product_stiffness_vector(neq,u_g,kp)
  endif
  ! assemble kp across the regions
  r = f-kp
  z = dprecon_g*r

  call scatter_and_assemble(neq,z,z_g)

  p = z

  ! pcg iteration
  pcg: do iter = 1,CG_MAXITER
    call scatter_and_assemble(neq,p,p_g)

    call product_stiffness_vector(neq,p_g,kp)

    rz = dot_product_all_proc(r,z_g)
    pkp = dot_product_all_proc(p_g,kp)

    ! avoids division by zero
    if (pkp /= zero) then
      alpha = rz / pkp
    else
      alpha = zero
    endif
    u_g = u_g + alpha * p_g

    maxp = maxvec(abs(p_g)); maxu = maxvec(abs(u_g))

    ! avoids division by zero
    if (maxu /= zero) then
      tol = abs(alpha) * maxp / maxu
    else
      tol = abs(alpha) * maxp
    endif

    !if (abs(alpha)*maxvec(abs(p_g))/maxvec(abs(u_g)) <= CG_TOL) then
    if (tol <= CG_TOL) then
      ! returns number of iterations
      if (present(niter)) niter = iter
      ! all done
      return
    endif

    r0 = r
    r = r - alpha * kp
    z = dprecon_g * r

    call scatter_and_assemble(neq,z,z_g)

    !beta = dot_product_all_proc(r,z_g)/rz ! Fletcher–Reeves
    beta = dot_product_all_proc(r-r0,z_g)/rz !  Polak–Ribière

    p = z + beta * p
  enddo pcg

  if (myid == 0) print *,'ERROR: PCG solver doesn''t converge! Tolerance:',tol, &
                         'alpha/maxp/maxu',alpha,maxp,maxu
  call synchronize_all()

  ! abort
  call exit_MPI(myrank,'PCG solver does not converge in routine diagpcg_solver()')

  end subroutine diagpcg_solver

!
!============================================
!

! diagonally preconditioned conjugate-gradient solver

  subroutine diagpcg_solver3(myid,neq,u_g,f,dprecon_g,niter)

  implicit none
  integer,intent(in) :: myid,neq
  real(kind=CUSTOM_REAL),dimension(0:neq),intent(inout) :: u_g
  real(kind=CUSTOM_REAL),dimension(0:neq),intent(in) :: f,dprecon_g
  integer, optional :: niter

  ! local parameters
  integer :: iter
  real(kind=CUSTOM_REAL) :: alpha,beta,pkp,rz,maxf,maxp,maxu,tol
  real(kind=CUSTOM_REAL),dimension(0:neq) :: kp,p,p_g,r0,r,z,z_g!,r_g

  real(kind=CUSTOM_REAL),parameter :: zero=0.0_CUSTOM_REAL,zerotol=1.0e-12_CUSTOM_REAL

  ! all global array variables are both MPI and regionally assembled.
  ! local array variables are regionally assembled.
  ! for MPI assembly of such array, we have to scatter again to region with
  ! regionally assembled values.

  ! PCG solver
  maxf = maxvec(abs(f))

  !debug
  !if (myid == 0) print *,'PCG solver: max load: ',maxf

  ! check if RHS is 0
  if (maxf <= zerotol) then
    u_g(:) = zero
    ! returns number of iterations
    if (present(niter)) niter = 0
    ! all done
    return
  endif

  kp = zero
  if (maxval(abs(u_g)) > zero) then
    call product_stiffness_vector3(neq,u_g,kp)
  endif
  ! assemble kp across the regions
  r = f-kp
  z = dprecon_g*r

  call scatter_and_assemble3(neq,z,z_g)

  p = z
  ! pcg iteration
  pcg: do iter = 1,CG_MAXITER
    call scatter_and_assemble3(neq,p,p_g)
    !call assemble_ghosts(myid,ngpart,maxngnode,nndof,neq,p,p_g) !,gdof)
    !print *,'pcg_bp4'

    call product_stiffness_vector3(neq,p_g,kp)

    rz = dot_product_all_proc(r,z_g)
    pkp = dot_product_all_proc(p_g,kp)

    ! avoids division by zero
    if (pkp /= zero) then
      alpha = rz / pkp
    else
      alpha = zero
    endif
    u_g = u_g + alpha * p_g

    maxp = maxvec(abs(p_g)); maxu = maxvec(abs(u_g))

    ! avoids division by zero
    if (maxu /= zero) then
      tol = abs(alpha) * maxp / maxu
    else
      tol = abs(alpha) * maxp
    endif

    if (tol <= CG_TOL) then
      ! returns number of iterations
      if (present(niter)) niter = iter
      ! all done
      return
    endif

    r0 = r !PR: Polak–Ribière
    r = r - alpha * kp
    z = dprecon_g * r

    !print *,'pcg_bp8'
    call scatter_and_assemble3(neq,z,z_g)
    beta = dot_product_all_proc(r-r0,z_g)/rz !PR: Polak–Ribière
    p = z + beta * p
  enddo pcg

  if (myid == 0) print *,'ERROR: PCG solver doesn''t converge! Tolerance:',tol, &
                         'alpha/maxp/maxu',alpha,maxp,maxu
  call synchronize_all()

  ! abort
  call exit_MPI(myrank,'PCG solver does not converge in routine diagpcg_solver3()')

  end subroutine diagpcg_solver3

!
!============================================
!

!TODO: this subroutine can be made non-blocking decomposing inner and outer elements

  subroutine product_stiffness_vector(neq,p_g,kp)

  use specfem_par, only: ADD_TRINF,NGLLCUBE,NSPEC_INNER_CORE,NSPEC_OUTER_CORE,NSPEC_CRUST_MANTLE, &
    NSPEC_TRINFINITE,NSPEC_INFINITE,NGLOB_INNER_CORE,NGLOB_OUTER_CORE,NGLOB_CRUST_MANTLE, &
    NGLOB_TRINFINITE,NGLOB_INFINITE

  use specfem_par_full_gravity, only: &
    gdof_cm,inode_elmt_cm,storekmat_crust_mantle, &
    gdof_oc,inode_elmt_oc,storekmat_outer_core, &
    gdof_ic,inode_elmt_ic,storekmat_inner_core, &
    gdof_trinf,inode_elmt_trinf,storekmat_trinfinite, &
    gdof_inf,inode_elmt_inf,storekmat_infinite

  implicit none
  integer,intent(in) :: neq
  real(kind=CUSTOM_REAL),intent(in) :: p_g(0:neq)
  real(kind=CUSTOM_REAL),intent(out) :: kp(0:neq)

  real(kind=CUSTOM_REAL) :: km(NGLLCUBE,NGLLCUBE),km_trinf(NGLLCUBE,NGLLCUBE),km_inf(NGLLCUBE,NGLLCUBE)
  real(kind=CUSTOM_REAL) :: kp_ic(NGLOB_INNER_CORE),kp_oc(NGLOB_OUTER_CORE), &
                            kp_cm(NGLOB_CRUST_MANTLE),kp_trinf(NGLOB_TRINFINITE),kp_inf(NGLOB_INFINITE)

  real(kind=CUSTOM_REAL),parameter :: zero=0.0_CUSTOM_REAL
  integer :: i_elmt,inode(NGLLCUBE),igdof(NGLLCUBE),inode_trinf(NGLLCUBE),igdof_trinf(NGLLCUBE), &
             inode_inf(NGLLCUBE),igdof_inf(NGLLCUBE)

  ! inner core
  kp_ic = zero
  do i_elmt = 1,NSPEC_INNER_CORE
    inode = inode_elmt_ic(:,i_elmt)
    igdof = gdof_ic(inode)
    km = storekmat_inner_core(:,:,i_elmt)
    kp_ic(inode) = kp_ic(inode)+matmul(km,p_g(igdof))
  enddo

  ! outer core
  kp_oc = zero
  do i_elmt = 1,NSPEC_OUTER_CORE
    inode = inode_elmt_oc(:,i_elmt)
    igdof = gdof_oc(inode)
    km = storekmat_outer_core(:,:,i_elmt)
    kp_oc(inode) = kp_oc(inode)+matmul(km,p_g(igdof))
  enddo

  ! crust mantle
  kp_cm = zero
  do i_elmt = 1,NSPEC_CRUST_MANTLE
    inode = inode_elmt_cm(:,i_elmt)
    igdof = gdof_cm(inode)
    km = storekmat_crust_mantle(:,:,i_elmt)
    kp_cm(inode) = kp_cm(inode)+matmul(km,p_g(igdof))
  enddo

  ! transition infinite
  if (ADD_TRINF) then
    kp_trinf = zero
    do i_elmt = 1,NSPEC_TRINFINITE
      inode_trinf = inode_elmt_trinf(:,i_elmt)
      igdof_trinf = gdof_trinf(inode_trinf)
      km_trinf = storekmat_trinfinite(:,:,i_elmt)
      kp_trinf(inode_trinf) = kp_trinf(inode_trinf)+matmul(km_trinf,p_g(igdof_trinf))
    enddo
  endif

  ! infinite
  kp_inf = zero
  do i_elmt = 1,NSPEC_INFINITE
    inode_inf = inode_elmt_inf(:,i_elmt)
    igdof_inf = gdof_inf(inode_inf)
    km_inf = storekmat_infinite(:,:,i_elmt)
    kp_inf(inode_inf) = kp_inf(inode_inf)+matmul(km_inf,p_g(igdof_inf))
  enddo

  ! assemble acroos the regions but not across the MPIs
  kp = zero

  ! crust_mantle
  kp(gdof_cm) = kp(gdof_cm)+kp_cm
  ! outer core
  kp(gdof_oc) = kp(gdof_oc)+kp_oc
  ! inner core
  kp(gdof_ic) = kp(gdof_ic)+kp_ic
  ! transition infinite
  if (ADD_TRINF) kp(gdof_trinf) = kp(gdof_trinf)+kp_trinf
  ! infinite
  kp(gdof_inf) = kp(gdof_inf)+kp_inf

  kp(0) = zero

  end subroutine product_stiffness_vector

!
!============================================
!

  subroutine scatter_and_assemble(neq,array,array_g)

  use specfem_par, only: ADD_TRINF,NPROCTOT_VAL,NGLOB_INNER_CORE,NGLOB_OUTER_CORE, &
    NGLOB_CRUST_MANTLE,NGLOB_TRINFINITE,NGLOB_INFINITE,NGLOB_CRUST_MANTLE

  use specfem_par, only: num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
    nibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle, &
    my_neighbors_crust_mantle, &
    num_interfaces_outer_core,max_nibool_interfaces_oc, &
    nibool_interfaces_outer_core,ibool_interfaces_outer_core, &
    my_neighbors_outer_core, &
    num_interfaces_inner_core,max_nibool_interfaces_ic, &
    nibool_interfaces_inner_core,ibool_interfaces_inner_core, &
    my_neighbors_inner_core

  use specfem_par_full_gravity, only: &
    num_interfaces_trinfinite,max_nibool_interfaces_trinfinite, &
    nibool_interfaces_trinfinite,ibool_interfaces_trinfinite,my_neighbors_trinfinite, &
    num_interfaces_infinite,max_nibool_interfaces_infinite, &
    nibool_interfaces_infinite,ibool_interfaces_infinite,my_neighbors_infinite

  use specfem_par_full_gravity, only: gdof_cm,gdof_oc,gdof_ic,gdof_trinf,gdof_inf

  implicit none
  integer,intent(in) :: neq
  real(kind=CUSTOM_REAL),intent(in) :: array(0:neq)
  real(kind=CUSTOM_REAL),intent(out) :: array_g(0:neq)

  real(kind=CUSTOM_REAL) :: array_ic(NGLOB_INNER_CORE),array_oc(NGLOB_OUTER_CORE), &
                            array_cm(NGLOB_CRUST_MANTLE),array_trinf(NGLOB_TRINFINITE),array_inf(NGLOB_INFINITE)

  real(kind=CUSTOM_REAL),parameter :: zero = 0.0_CUSTOM_REAL

  ! scatter array
  array_ic = array(gdof_ic)
  array_oc = array(gdof_oc)
  array_cm = array(gdof_cm)
  if (ADD_TRINF) array_trinf = array(gdof_trinf)
  array_inf = array(gdof_inf)

  ! assemble across the MPI processes in a region
  ! crust_mantle
  call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_CRUST_MANTLE,array_cm, &
                           num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                           nibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle, &
                           my_neighbors_crust_mantle)

  ! outer core
  call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_OUTER_CORE,array_oc, &
                           num_interfaces_outer_core,max_nibool_interfaces_oc, &
                           nibool_interfaces_outer_core,ibool_interfaces_outer_core, &
                           my_neighbors_outer_core)

  ! inner core
  call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_INNER_CORE,array_ic, &
                           num_interfaces_inner_core,max_nibool_interfaces_ic, &
                           nibool_interfaces_inner_core,ibool_interfaces_inner_core, &
                           my_neighbors_inner_core)

  ! transition infinite
  if (ADD_TRINF) then
    call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_TRINFINITE,array_trinf, &
                             num_interfaces_trinfinite,max_nibool_interfaces_trinfinite, &
                             nibool_interfaces_trinfinite,ibool_interfaces_trinfinite,my_neighbors_trinfinite)
  endif

  ! infinite
  call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_INFINITE,array_inf, &
                           num_interfaces_infinite,max_nibool_interfaces_infinite, &
                           nibool_interfaces_infinite,ibool_interfaces_infinite,my_neighbors_infinite)

  ! gather from all regions but not assemble since it is already assembled across
  array_g = zero

  ! the regions assemble across the different regions in a process
  ! crust_mantle
  array_g(gdof_cm) = array_cm
  ! outer core
  array_g(gdof_oc) = array_oc
  ! inner core
  array_g(gdof_ic) = array_ic
  ! transition infinite
  if (ADD_TRINF) array_g(gdof_trinf) = array_trinf
  ! infinite
  array_g(gdof_inf) = array_inf

  array_g(0) = zero

  end subroutine scatter_and_assemble

!
!============================================
!

! TODO: this subroutine can be made non-blocking decomposing inner and outer elements

  subroutine product_stiffness_vector3(neq,p_g,kp)

  use constants, only: ADD_TRINF,NGLLCUBE_INF,IFLAG_IN_FICTITIOUS_CUBE

  use specfem_par, only: NSPEC_INNER_CORE,NSPEC_OUTER_CORE,NSPEC_CRUST_MANTLE, &
    NSPEC_TRINFINITE,NSPEC_INFINITE

  !use specfem_par_innercore, only: idoubling_inner_core

  use specfem_par_full_gravity, only: &
    nnode_ic1,nnode_oc1,nnode_cm1,nnode_trinf1,nnode_inf1, &
    gdof_cm1,inode_elmt_cm1,storekmat_crust_mantle1, &
    gdof_oc1,inode_elmt_oc1,storekmat_outer_core1, &
    gdof_ic1,inode_elmt_ic1,storekmat_inner_core1, &
    gdof_trinf1,inode_elmt_trinf1,storekmat_trinfinite1, &
    gdof_inf1,inode_elmt_inf1,storekmat_infinite1

  implicit none
  integer,intent(in) :: neq
  real(kind=CUSTOM_REAL),intent(in) :: p_g(0:neq)
  real(kind=CUSTOM_REAL),intent(out) :: kp(0:neq)

  real(kind=CUSTOM_REAL) :: km(NGLLCUBE_INF,NGLLCUBE_INF),km_trinf(NGLLCUBE_INF,NGLLCUBE_INF),km_inf(NGLLCUBE_INF,NGLLCUBE_INF)
  real(kind=CUSTOM_REAL) :: kp_ic(nnode_ic1),kp_oc(nnode_oc1),kp_cm(nnode_cm1), &
                            kp_trinf(nnode_trinf1),kp_inf(nnode_inf1)

  real(kind=CUSTOM_REAL),parameter :: zero=0.0_CUSTOM_REAL
  integer :: i_elmt
  integer :: inode(NGLLCUBE_INF),igdof(NGLLCUBE_INF),inode_trinf(NGLLCUBE_INF),igdof_trinf(NGLLCUBE_INF), &
             inode_inf(NGLLCUBE_INF),igdof_inf(NGLLCUBE_INF)

  ! explicit loop
  !integer :: i,j,iglob,igll
  !real(kind=CUSTOM_REAL) :: sum
  !real(kind=CUSTOM_REAL) :: p_g_elmt(NGLLCUBE_INF)

  ! inner core
  kp_ic(:) = zero
  do i_elmt = 1,NSPEC_INNER_CORE
    !? if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE) cycle
    inode(:) = inode_elmt_ic1(:,i_elmt)
    igdof(:) = gdof_ic1(inode(:))
    km = storekmat_inner_core1(:,:,i_elmt)
    kp_ic(inode(:)) = kp_ic(inode(:)) + matmul(km,p_g(igdof))
  enddo

  ! outer core
  kp_oc(:) = zero
  do i_elmt = 1,NSPEC_OUTER_CORE
    inode(:) = inode_elmt_oc1(:,i_elmt)
    igdof(:) = gdof_oc1(inode(:))
    km = storekmat_outer_core1(:,:,i_elmt)
    kp_oc(inode(:)) = kp_oc(inode) + matmul(km,p_g(igdof))
  enddo

  ! crust mantle
  kp_cm(:) = zero
  do i_elmt = 1,NSPEC_CRUST_MANTLE
    inode(:) = inode_elmt_cm1(:,i_elmt)
    igdof(:) = gdof_cm1(inode(:))

    ! w/ intrinsic matmul: background solve -> 11.40 s
    km = storekmat_crust_mantle1(:,:,i_elmt)
    kp_cm(inode(:)) = kp_cm(inode(:)) + matmul(km,p_g(igdof))

    ! w/out km copy: background solve -> 13.74 s
    !kp_cm(inode(:)) = kp_cm(inode(:)) + matmul(storekmat_crust_mantle1(:,:,i_elmt),p_g(igdof))

    ! w/ local copies: background solve -> 11.48 s
    !km = storekmat_crust_mantle1(:,:,i_elmt)
    !p_g_elmt(:) = p_g(igdof(:))
    !kp_cm(inode(:)) = kp_cm(inode(:)) + matmul(km,p_g_elmt)

    ! w/ explicit loop: background solve -> 14.88 s
    ! with kp_i = km_ij * pg_j    with kp = (n), kmat = (n x n), pg = (n)
    !do i = 1,NGLLCUBE_INF
    !  sum = 0.0_CUSTOM_REAL
    !  do j = 1,NGLLCUBE_INF
    !    iglob = igdof(j)
    !    sum = sum + km(i,j) * p_g(iglob)
    !  enddo
    !  igll = inode(i)
    !  kp_cm(igll) = kp_cm(igll) + sum
    !enddo
  enddo

  ! transition infinite
  if (ADD_TRINF) then
    kp_trinf(:) = zero
    do i_elmt = 1,NSPEC_TRINFINITE
      inode_trinf(:) = inode_elmt_trinf1(:,i_elmt)
      igdof_trinf(:) = gdof_trinf1(inode_trinf(:))
      km_trinf = storekmat_trinfinite1(:,:,i_elmt)
      kp_trinf(inode_trinf(:)) = kp_trinf(inode_trinf(:)) + matmul(km_trinf,p_g(igdof_trinf))
    enddo
  endif

  ! infinite
  kp_inf(:) = zero
  do i_elmt = 1,NSPEC_INFINITE
    inode_inf(:) = inode_elmt_inf1(:,i_elmt)
    igdof_inf(:) = gdof_inf1(inode_inf(:))
    km_inf = storekmat_infinite1(:,:,i_elmt)
    kp_inf(inode_inf(:)) = kp_inf(inode_inf(:)) + matmul(km_inf,p_g(igdof_inf))
  enddo

  ! assemble acroos the regions but not across the MPIs
  kp(:) = zero

  ! assemble across the different regions in a process
  ! crust_mantle
  kp(gdof_cm1(:)) = kp(gdof_cm1(:)) + kp_cm(:)
  ! outer core
  kp(gdof_oc1(:)) = kp(gdof_oc1(:)) + kp_oc(:)
  ! inner core
  kp(gdof_ic1(:)) = kp(gdof_ic1(:)) + kp_ic(:)
  ! transition infinite
  if (ADD_TRINF) kp(gdof_trinf1(:)) = kp(gdof_trinf1(:)) + kp_trinf(:)
  ! infinite
  kp(gdof_inf1(:)) = kp(gdof_inf1(:)) + kp_inf(:)

  kp(0) = zero

  end subroutine product_stiffness_vector3

!
!============================================
!

  subroutine scatter_and_assemble3(neq,array,array_g)

  use specfem_par, only: ADD_TRINF,NPROCTOT_VAL

  use specfem_par_full_gravity, only: num_interfaces_crust_mantle1,max_nibool_interfaces_crust_mantle1, &
    nibool_interfaces_crust_mantle1,ibool_interfaces_crust_mantle1, &
    my_neighbors_crust_mantle1, &
    num_interfaces_outer_core1,max_nibool_interfaces_outer_core1, &
    nibool_interfaces_outer_core1,ibool_interfaces_outer_core1, &
    my_neighbors_outer_core1, &
    num_interfaces_inner_core1,max_nibool_interfaces_inner_core1, &
    nibool_interfaces_inner_core1,ibool_interfaces_inner_core1, &
    my_neighbors_inner_core1, &
    num_interfaces_trinfinite1,max_nibool_interfaces_trinfinite1, &
    nibool_interfaces_trinfinite1,ibool_interfaces_trinfinite1, &
    my_neighbors_trinfinite1, &
    num_interfaces_infinite1,max_nibool_interfaces_infinite1, &
    nibool_interfaces_infinite1,ibool_interfaces_infinite1, &
    my_neighbors_infinite1, &
    nnode_ic1,nnode_oc1,nnode_cm1,nnode_trinf1,nnode_inf1

  use specfem_par_full_gravity, only: gdof_cm1,gdof_oc1,gdof_ic1,gdof_trinf1,gdof_inf1

  implicit none
  integer,intent(in) :: neq
  real(kind=CUSTOM_REAL),intent(in) :: array(0:neq)
  real(kind=CUSTOM_REAL),intent(out) :: array_g(0:neq)

  ! local parameters
  real(kind=CUSTOM_REAL) :: array_ic(nnode_ic1),array_oc(nnode_oc1),array_cm(nnode_cm1), &
                            array_trinf(nnode_trinf1),array_inf(nnode_inf1)

  real(kind=CUSTOM_REAL),parameter :: zero = 0.0_CUSTOM_REAL

  ! scatter array
  array_ic(:) = array(gdof_ic1(:))
  array_oc(:) = array(gdof_oc1(:))
  array_cm(:) = array(gdof_cm1(:))
  if (ADD_TRINF) array_trinf(:) = array(gdof_trinf1(:))
  array_inf(:) = array(gdof_inf1(:))

  ! assemble across the MPI processes in a region
  ! crust_mantle
  call assemble_MPI_scalar(NPROCTOT_VAL,nnode_cm1,array_cm, &
                           num_interfaces_crust_mantle1,max_nibool_interfaces_crust_mantle1, &
                           nibool_interfaces_crust_mantle1,ibool_interfaces_crust_mantle1, &
                           my_neighbors_crust_mantle1)

  ! outer core
  call assemble_MPI_scalar(NPROCTOT_VAL,nnode_oc1,array_oc, &
                           num_interfaces_outer_core1,max_nibool_interfaces_outer_core1, &
                           nibool_interfaces_outer_core1,ibool_interfaces_outer_core1, &
                           my_neighbors_outer_core1)

  ! inner core
  call assemble_MPI_scalar(NPROCTOT_VAL,nnode_ic1,array_ic, &
                           num_interfaces_inner_core1,max_nibool_interfaces_inner_core1, &
                           nibool_interfaces_inner_core1,ibool_interfaces_inner_core1, &
                           my_neighbors_inner_core1)

  ! transition infinite
  if (ADD_TRINF) then
    call assemble_MPI_scalar(NPROCTOT_VAL,nnode_trinf1,array_trinf, &
                             num_interfaces_trinfinite1,max_nibool_interfaces_trinfinite1, &
                             nibool_interfaces_trinfinite1,ibool_interfaces_trinfinite1,my_neighbors_trinfinite1)
  endif

  ! infinite
  call assemble_MPI_scalar(NPROCTOT_VAL,nnode_inf1,array_inf, &
                           num_interfaces_infinite1,max_nibool_interfaces_infinite1, &
                           nibool_interfaces_infinite1,ibool_interfaces_infinite1,my_neighbors_infinite1)

  ! gather from all regions but not assemble since it is already assembled across
  ! the regions assemble across the different regions in a process
  array_g(:) = zero

  ! crust_mantle
  array_g(gdof_cm1(:)) = array_cm(:)
  ! outer core
  array_g(gdof_oc1(:)) = array_oc(:)
  ! inner core
  array_g(gdof_ic1(:)) = array_ic(:)
  ! transition infinite
  if (ADD_TRINF) array_g(gdof_trinf1(:)) = array_trinf(:)
  ! infinite
  array_g(gdof_inf1(:)) = array_inf(:)

  array_g(0) = zero

  end subroutine scatter_and_assemble3

!
!============================================
!

! interpolate solution to original mesh to compute initial guess

  subroutine interpolate3to5(nelmt,nnode,nnode1,inode_elmt,nmir,inode_map,is_active_gll,igll_active_on,x3,x5)

  use constants, only: NDIM,NGLLX,NGLLX_INF,NGLLCUBE,NGLLCUBE_INF

  use siem_gll_library, only: gll_quadrature3inNGLL,zwgljd

  implicit none
  integer,intent(in) :: nelmt,nnode,nnode1
  integer,intent(in) :: inode_elmt(NGLLCUBE,nelmt)
  integer,intent(in) :: nmir(nnode)
  integer,intent(in) :: inode_map(2,nnode)
  integer,intent(in) :: igll_active_on(NGLLCUBE_INF)
  logical,intent(in) :: is_active_gll(NGLLCUBE)

  real(kind=CUSTOM_REAL),intent(in) :: x3(nnode1) ! array for 3 GLLX points
  real(kind=CUSTOM_REAL),intent(out) :: x5(nnode) ! array for 5 GLLX points

  ! local parameters
  double precision :: lagrange_gll3inNGLL(NGLLCUBE,27)
  double precision :: xigll(NGLLX),wxgll(NGLLX)
  double precision :: xigll1(NGLLX_INF),wxgll1(NGLLX_INF)

  integer :: i_node,ielmt,igll,inode1
  integer :: inodes1(NGLLCUBE_INF)

  call zwgljd(xigll1,wxgll1,NGLLX_INF,0.d0,0.d0)
  call zwgljd(xigll,wxgll,NGLLX,0.d0,0.d0)

  call gll_quadrature3inNGLL(NDIM,NGLLX,NGLLCUBE,xigll,xigll1,lagrange_gll3inNGLL)

  x5(:) = 0.0_CUSTOM_REAL

  do i_node = 1,nnode
    ielmt = inode_map(1,i_node)

    ! skip fictitious nodes
    if (ielmt <= 0) cycle

    igll = inode_map(2,i_node)
    if (is_active_gll(igll)) then
      ! takes value of mirrored node in Level-1 array
      inode1 = nmir(i_node)
      x5(i_node) = x3(inode1)
    else
      ! interpolate values
      inodes1(:) = nmir(inode_elmt(igll_active_on(:),ielmt))
      x5(i_node) = real(sum(lagrange_gll3inNGLL(igll,:)*x3(inodes1(:))),kind=CUSTOM_REAL)
    endif
  enddo

  end subroutine interpolate3to5


end module siem_solver_mpi

