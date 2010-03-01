!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            March 2010
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

 
  subroutine save_kernels_crust_mantle(myrank,scale_t,scale_displ, &
                  cijkl_kl_crust_mantle,rho_kl_crust_mantle, &
                  alpha_kl_crust_mantle,beta_kl_crust_mantle, &  
                  ystore_crust_mantle,zstore_crust_mantle, &
                  rhostore_crust_mantle,muvstore_crust_mantle, &
                  kappavstore_crust_mantle,ibool_crust_mantle, &
                  LOCAL_PATH)

  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer myrank
  
  double precision :: scale_t,scale_displ

  real(kind=CUSTOM_REAL), dimension(21,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
    cijkl_kl_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
    rho_kl_crust_mantle, beta_kl_crust_mantle, alpha_kl_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: &
        ystore_crust_mantle,zstore_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ISO_MANTLE) :: &
        rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle
    
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle

  character(len=150) LOCAL_PATH

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
    mu_kl_crust_mantle, kappa_kl_crust_mantle, rhonotprime_kl_crust_mantle
  real(kind=CUSTOM_REAL),dimension(21) ::  cijkl_kl_local   
  real(kind=CUSTOM_REAL) :: scale_kl,scale_kl_ani,scale_kl_rho
  real(kind=CUSTOM_REAL) :: rhol,mul,kappal,rho_kl,alpha_kl,beta_kl
  integer :: ispec,i,j,k,iglob
  character(len=150) prname
  

  scale_kl = scale_t/scale_displ * 1.d9
  ! For anisotropic kernels
  ! final unit : [s km^(-3) GPa^(-1)]
  scale_kl_ani = scale_t**3 / (RHOAV*R_EARTH**3) * 1.d18
  ! final unit : [s km^(-3) (kg/m^3)^(-1)]
  scale_kl_rho = scale_t / scale_displ / RHOAV * 1.d9

  ! crust_mantle
  do ispec = 1, NSPEC_CRUST_MANTLE
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX


          if (ANISOTROPIC_KL) then

            ! For anisotropic kernels
            iglob = ibool_crust_mantle(i,j,k,ispec)

            ! The cartesian global cijkl_kl are rotated into the spherical local cijkl_kl
            ! ystore and zstore are thetaval and phival (line 2252) -- dangerous
            call rotate_kernels_dble(cijkl_kl_crust_mantle(:,i,j,k,ispec),cijkl_kl_local, &
                 ystore_crust_mantle(iglob),zstore_crust_mantle(iglob))
                 
            cijkl_kl_crust_mantle(:,i,j,k,ispec) = cijkl_kl_local * scale_kl_ani
            rho_kl_crust_mantle(i,j,k,ispec) = rho_kl_crust_mantle(i,j,k,ispec) * scale_kl_rho

          else

            rhol = rhostore_crust_mantle(i,j,k,ispec)
            mul = muvstore_crust_mantle(i,j,k,ispec)
            kappal = kappavstore_crust_mantle(i,j,k,ispec)
            
            ! kernel values for rho, kappa, mu
            rho_kl = - rhol * rho_kl_crust_mantle(i,j,k,ispec)
            alpha_kl = - kappal * alpha_kl_crust_mantle(i,j,k,ispec)
            beta_kl =  - 2 * mul * beta_kl_crust_mantle(i,j,k,ispec)

            rhonotprime_kl_crust_mantle(i,j,k,ispec) = rho_kl * scale_kl
            mu_kl_crust_mantle(i,j,k,ispec) = beta_kl * scale_kl
            kappa_kl_crust_mantle(i,j,k,ispec) = alpha_kl * scale_kl
            
            ! kernels rho^prime, beta, alpha
            rho_kl_crust_mantle(i,j,k,ispec) = (rho_kl + alpha_kl + beta_kl) * scale_kl
            beta_kl_crust_mantle(i,j,k,ispec) = 2 * (beta_kl - FOUR_THIRDS * mul * alpha_kl / kappal) * scale_kl
            alpha_kl_crust_mantle(i,j,k,ispec) = 2 * (1 +  FOUR_THIRDS * mul / kappal) * alpha_kl * scale_kl

          endif

        enddo
      enddo
    enddo
  enddo

  call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_PATH)

  ! For anisotropic kernels
  if (ANISOTROPIC_KL) then

    open(unit=27,file=trim(prname)//'rho_kernel.bin',status='unknown',form='unformatted',action='write')
    write(27) -rho_kl_crust_mantle
    close(27)
    open(unit=27,file=trim(prname)//'cijkl_kernel.bin',status='unknown',form='unformatted',action='write')
    write(27) -cijkl_kl_crust_mantle
    close(27)

  else

    open(unit=27,file=trim(prname)//'rhonotprime_kernel.bin',status='unknown',form='unformatted',action='write')
    write(27) rhonotprime_kl_crust_mantle
    close(27)
    open(unit=27,file=trim(prname)//'kappa_kernel.bin',status='unknown',form='unformatted',action='write')
    write(27) kappa_kl_crust_mantle
    close(27)
    open(unit=27,file=trim(prname)//'mu_kernel.bin',status='unknown',form='unformatted',action='write')
    write(27) mu_kl_crust_mantle
    close(27)


    open(unit=27,file=trim(prname)//'rho_kernel.bin',status='unknown',form='unformatted',action='write')
    write(27) rho_kl_crust_mantle
    close(27)
    open(unit=27,file=trim(prname)//'alpha_kernel.bin',status='unknown',form='unformatted',action='write')
    write(27) alpha_kl_crust_mantle
    close(27)
    open(unit=27,file=trim(prname)//'beta_kernel.bin',status='unknown',form='unformatted',action='write')
    write(27) beta_kl_crust_mantle
    close(27)

  endif

 
  end subroutine save_kernels_crust_mantle

!  
!-------------------------------------------------------------------------------------------------
!  
 
  subroutine save_kernels_outer_core(myrank,scale_t,scale_displ, &
                        rho_kl_outer_core,alpha_kl_outer_core, &
                        rhostore_outer_core,kappavstore_outer_core, &
                        deviatoric_outercore,nspec_beta_kl_outer_core,beta_kl_outer_core, &
                        LOCAL_PATH)
                        
  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer myrank
  
  double precision :: scale_t,scale_displ

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ADJOINT) :: &
    rho_kl_outer_core,alpha_kl_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: &
        rhostore_outer_core,kappavstore_outer_core

  integer nspec_beta_kl_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_beta_kl_outer_core) :: &
    beta_kl_outer_core
  logical deviatoric_outercore

  character(len=150) LOCAL_PATH

  ! local parameters
  real(kind=CUSTOM_REAL):: scale_kl
  real(kind=CUSTOM_REAL) :: rhol,kappal,rho_kl,alpha_kl,beta_kl  
  integer :: ispec,i,j,k
  character(len=150) prname
  
  scale_kl = scale_t/scale_displ * 1.d9

  ! outer_core
  do ispec = 1, NSPEC_OUTER_CORE
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          rhol = rhostore_outer_core(i,j,k,ispec)
          kappal = kappavstore_outer_core(i,j,k,ispec)
          rho_kl = - rhol * rho_kl_outer_core(i,j,k,ispec)
          alpha_kl = - kappal * alpha_kl_outer_core(i,j,k,ispec)
          
          rho_kl_outer_core(i,j,k,ispec) = (rho_kl + alpha_kl) * scale_kl
          alpha_kl_outer_core(i,j,k,ispec) = 2 * alpha_kl * scale_kl
          
          
          !deviatoric kernel check
          if( deviatoric_outercore ) then
            beta_kl =  - 2 * beta_kl_outer_core(i,j,k,ispec)  ! not using mul, since it's zero for the fluid
            beta_kl_outer_core(i,j,k,ispec) = beta_kl
          endif
          
        enddo
      enddo
    enddo
  enddo

  call create_name_database(prname,myrank,IREGION_OUTER_CORE,LOCAL_PATH)
  open(unit=27,file=trim(prname)//'rho_kernel.bin',status='unknown',form='unformatted',action='write')
  write(27) rho_kl_outer_core
  close(27)
  open(unit=27,file=trim(prname)//'alpha_kernel.bin',status='unknown',form='unformatted',action='write')
  write(27) alpha_kl_outer_core
  close(27)

  !deviatoric kernel check
  if( deviatoric_outercore ) then
    open(unit=27,file=trim(prname)//'mu_kernel.bin',status='unknown',form='unformatted',action='write')
    write(27) beta_kl_outer_core
    close(27)        
  endif
 
  end subroutine save_kernels_outer_core
  
!  
!-------------------------------------------------------------------------------------------------
!  
 
  subroutine save_kernels_inner_core(myrank,scale_t,scale_displ, &
                          rho_kl_inner_core,beta_kl_inner_core,alpha_kl_inner_core, &
                          rhostore_inner_core,muvstore_inner_core,kappavstore_inner_core, &
                          LOCAL_PATH)
  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer myrank
  
  double precision :: scale_t,scale_displ
  
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ADJOINT) :: &
    rho_kl_inner_core, beta_kl_inner_core, alpha_kl_inner_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: &
        rhostore_inner_core, kappavstore_inner_core,muvstore_inner_core

  character(len=150) LOCAL_PATH

  ! local parameters
  real(kind=CUSTOM_REAL):: scale_kl
  real(kind=CUSTOM_REAL) :: rhol,mul,kappal,rho_kl,alpha_kl,beta_kl  
  integer :: ispec,i,j,k
  character(len=150) prname
  

  scale_kl = scale_t/scale_displ * 1.d9

  ! inner_core
  do ispec = 1, NSPEC_INNER_CORE
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          rhol = rhostore_inner_core(i,j,k,ispec)
          mul = muvstore_inner_core(i,j,k,ispec)
          kappal = kappavstore_inner_core(i,j,k,ispec)
          
          rho_kl = -rhol * rho_kl_inner_core(i,j,k,ispec)
          alpha_kl = -kappal * alpha_kl_inner_core(i,j,k,ispec)
          beta_kl =  - 2 * mul * beta_kl_inner_core(i,j,k,ispec)
          
          rho_kl_inner_core(i,j,k,ispec) = (rho_kl + alpha_kl + beta_kl) * scale_kl
          beta_kl_inner_core(i,j,k,ispec) = 2 * (beta_kl - FOUR_THIRDS * mul * alpha_kl / kappal) * scale_kl
          alpha_kl_inner_core(i,j,k,ispec) = 2 * (1 +  FOUR_THIRDS * mul / kappal) * alpha_kl * scale_kl
        enddo
      enddo
    enddo
  enddo

  call create_name_database(prname,myrank,IREGION_INNER_CORE,LOCAL_PATH)
  open(unit=27,file=trim(prname)//'rho_kernel.bin',status='unknown',form='unformatted',action='write')
  write(27) rho_kl_inner_core
  close(27)
  open(unit=27,file=trim(prname)//'alpha_kernel.bin',status='unknown',form='unformatted',action='write')
  write(27) alpha_kl_inner_core
  close(27)
  open(unit=27,file=trim(prname)//'beta_kernel.bin',status='unknown',form='unformatted',action='write')
  write(27) beta_kl_inner_core
  close(27)
  
  end subroutine save_kernels_inner_core

!  
!-------------------------------------------------------------------------------------------------
!  
 
  subroutine save_kernels_boundary_kl(myrank,scale_t,scale_displ, &
                                  moho_kl,d400_kl,d670_kl,cmb_kl,icb_kl, &
                                  LOCAL_PATH,HONOR_1D_SPHERICAL_MOHO)

  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer myrank
  
  double precision :: scale_t,scale_displ

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_MOHO) :: moho_kl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_400) :: d400_kl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_670) :: d670_kl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_CMB) :: cmb_kl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_ICB) :: icb_kl

  character(len=150) LOCAL_PATH

  logical HONOR_1D_SPHERICAL_MOHO
  
  ! local parameters
  real(kind=CUSTOM_REAL):: scale_kl
  character(len=150) prname


  scale_kl = scale_t/scale_displ * 1.d9

  ! scale the boundary kernels properly: *scale_kl gives s/km^3 and 1.d3 gives
  ! the relative boundary kernels (for every 1 km) in s/km^2
  moho_kl = moho_kl * scale_kl * 1.d3
  d400_kl = d400_kl * scale_kl * 1.d3
  d670_kl = d670_kl * scale_kl * 1.d3
  cmb_kl = cmb_kl * scale_kl * 1.d3
  icb_kl = icb_kl * scale_kl * 1.d3

  call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_PATH)
  
  if (.not. SUPPRESS_CRUSTAL_MESH .and. HONOR_1D_SPHERICAL_MOHO) then
    open(unit=27,file=trim(prname)//'moho_kernel.bin',status='unknown',form='unformatted',action='write')
    write(27) moho_kl
    close(27)
  endif
  
  open(unit=27,file=trim(prname)//'d400_kernel.bin',status='unknown',form='unformatted',action='write')
  write(27) d400_kl
  close(27)
  
  open(unit=27,file=trim(prname)//'d670_kernel.bin',status='unknown',form='unformatted',action='write')
  write(27) d670_kl
  close(27)
  
  open(unit=27,file=trim(prname)//'CMB_kernel.bin',status='unknown',form='unformatted',action='write')
  write(27) cmb_kl
  close(27)
  
  call create_name_database(prname,myrank,IREGION_OUTER_CORE,LOCAL_PATH)
  
  open(unit=27,file=trim(prname)//'ICB_kernel.bin',status='unknown',form='unformatted',action='write')
  write(27) icb_kl
  close(27)

  
  end subroutine save_kernels_boundary_kl
  

!  
!-------------------------------------------------------------------------------------------------
!  
 
  subroutine save_kernels_source_derivatives(nrec_local,NSOURCES,scale_displ,scale_t, &
                                nu_source,moment_der,sloc_der,number_receiver_global)

  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer nrec_local,NSOURCES
  double precision :: scale_displ,scale_t

  double precision :: nu_source(NDIM,NDIM,NSOURCES)
  real(kind=CUSTOM_REAL) :: moment_der(NDIM,NDIM,nrec_local),sloc_der(NDIM,nrec_local)

  integer, dimension(nrec_local) :: number_receiver_global
  
  ! local parameters
  real(kind=CUSTOM_REAL),parameter :: scale_mass = RHOAV * (R_EARTH**3)
  integer :: irec_local
  character(len=150) outputname

  !scale_mass = RHOAV * (R_EARTH**3)

  do irec_local = 1, nrec_local
    ! rotate and scale the location derivatives to correspond to dn,de,dz
    sloc_der(:,irec_local) = matmul(nu_source(:,:,irec_local),sloc_der(:,irec_local)) &
                             * scale_displ * scale_t
                             
    ! rotate scale the moment derivatives to correspond to M[n,e,z][n,e,z]
    moment_der(:,:,irec_local) = matmul(matmul(nu_source(:,:,irec_local),moment_der(:,:,irec_local)),&
               transpose(nu_source(:,:,irec_local))) * scale_t ** 3 / scale_mass

    write(outputname,'(a,i5.5)') 'OUTPUT_FILES/src_frechet.',number_receiver_global(irec_local)
    open(unit=27,file=trim(outputname),status='unknown',action='write')
  !
  ! r -> z, theta -> -n, phi -> e, plus factor 2 for Mrt,Mrp,Mtp, and 1e-7 to dyne.cm
  !  Mrr =  Mzz
  !  Mtt =  Mnn
  !  Mpp =  Mee
  !  Mrt = -Mzn
  !  Mrp =  Mze
  !  Mtp = -Mne
  ! minus sign for sloc_der(3,irec_local) to get derivative for depth instead of radius

    write(27,'(g16.5)') moment_der(3,3,irec_local) * 1e-7
    write(27,'(g16.5)') moment_der(1,1,irec_local) * 1e-7
    write(27,'(g16.5)') moment_der(2,2,irec_local) * 1e-7
    write(27,'(g16.5)') -2*moment_der(1,3,irec_local) * 1e-7
    write(27,'(g16.5)') 2*moment_der(2,3,irec_local) * 1e-7
    write(27,'(g16.5)') -2*moment_der(1,2,irec_local) * 1e-7
    write(27,'(g16.5)') sloc_der(2,irec_local)
    write(27,'(g16.5)') sloc_der(1,irec_local)
    write(27,'(g16.5)') -sloc_der(3,irec_local)
    close(27)
  enddo


  end subroutine save_kernels_source_derivatives

    