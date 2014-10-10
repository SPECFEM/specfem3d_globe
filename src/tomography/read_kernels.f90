!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
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


subroutine read_kernels_iso()

! reads in smoothed kernels: bulk, beta, rho

  use tomography_kernels_iso

  implicit none
  real(kind=CUSTOM_REAL) :: min_vp,min_vs,max_vp,max_vs,min_rho,max_rho
  integer :: ier
  character(len=150) :: m_file, fname

  ! reads in smoothed (& summed) event kernel
  if (USE_ALPHA_BETA_RHO) then
    ! reads in alpha kernel
    fname = 'alpha_kernel_smooth'
  else
    ! reads in bulk_c kernel
    fname = 'bulk_c_kernel_smooth'
  endif
  write(m_file,'(a,i6.6,a)') 'INPUT_GRADIENT/proc',myrank,'_reg1_'//trim(fname)//'.bin'
  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print*,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_bulk(:,:,:,1:nspec)
  close(IIN)

  ! beta kernel
  if (USE_ALPHA_BETA_RHO) then
    ! reads in beta kernel
    fname = 'beta_kernel_smooth'
  else
    ! reads in bulk_beta kernel
    fname = 'bulk_beta_kernel_smooth'
  endif
  write(m_file,'(a,i6.6,a)') 'INPUT_GRADIENT/proc',myrank,'_reg1_'//trim(fname)//'.bin'
  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print*,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_beta(:,:,:,1:nspec)
  close(IIN)

  ! rho kernel
  if (USE_RHO_SCALING) then

    ! uses scaling relation with shear perturbations
    kernel_rho(:,:,:,:) = RHO_SCALING * kernel_beta(:,:,:,:)

  else

    ! uses rho kernel
    write(m_file,'(a,i6.6,a)') 'INPUT_GRADIENT/proc',myrank,'_reg1_rho_kernel_smooth.bin'
    open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print*,'Error opening: ',trim(m_file)
      call exit_mpi(myrank,'file not found')
    endif
    read(IIN) kernel_rho(:,:,:,1:nspec)
    close(IIN)
  endif

  ! statistics
  call min_all_cr(minval(kernel_bulk),min_vp)
  call max_all_cr(maxval(kernel_bulk),max_vp)

  call min_all_cr(minval(kernel_beta),min_vs)
  call max_all_cr(maxval(kernel_beta),max_vs)

  call min_all_cr(minval(kernel_rho),min_rho)
  call max_all_cr(maxval(kernel_rho),max_rho)

  if (myrank == 0) then
    print*,'initial kernels:'
    if (USE_ALPHA_BETA_RHO) then
      print*,'  alpha min/max: ',min_vp,max_vp
      print*,'  beta min/max : ',min_vs,max_vs
    else
      print*,'  bulk_c min/max   : ',min_vp,max_vp
      print*,'  bulk_beta min/max: ',min_vs,max_vs
    endif
    print*,'  rho min/max: ',min_rho,max_rho
    print*
  endif

end subroutine read_kernels_iso

!
!-------------------------------------------------------------------------------------------------
!


subroutine read_kernels_tiso()

! reads in smoothed kernels: bulk, betav, betah, eta

  use tomography_kernels_tiso

  implicit none
  real(kind=CUSTOM_REAL) :: min_vsv,min_vsh,max_vsv,max_vsh,min_eta,max_eta,min_bulk,max_bulk
  integer :: ier
  character(len=150) :: m_file

  ! bulk kernel
  write(m_file,'(a,i6.6,a)') 'INPUT_GRADIENT/proc',myrank,'_reg1_bulk_c_kernel_smooth.bin'
  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print*,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_bulk(:,:,:,1:nspec)
  close(IIN)

  ! betav kernel
  write(m_file,'(a,i6.6,a)') 'INPUT_GRADIENT/proc',myrank,'_reg1_bulk_betav_kernel_smooth.bin'
  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print*,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_betav(:,:,:,1:nspec)
  close(IIN)

  ! betah kernel
  write(m_file,'(a,i6.6,a)') 'INPUT_GRADIENT/proc',myrank,'_reg1_bulk_betah_kernel_smooth.bin'
  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print*,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_betah(:,:,:,1:nspec)
  close(IIN)

  ! eta kernel
  write(m_file,'(a,i6.6,a)') 'INPUT_GRADIENT/proc',myrank,'_reg1_eta_kernel_smooth.bin'
  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print*,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_eta(:,:,:,1:nspec)
  close(IIN)


  ! statistics
  call min_all_cr(minval(kernel_bulk),min_bulk)
  call max_all_cr(maxval(kernel_bulk),max_bulk)

  call min_all_cr(minval(kernel_betah),min_vsh)
  call max_all_cr(maxval(kernel_betah),max_vsh)

  call min_all_cr(minval(kernel_betav),min_vsv)
  call max_all_cr(maxval(kernel_betav),max_vsv)

  call min_all_cr(minval(kernel_eta),min_eta)
  call max_all_cr(maxval(kernel_eta),max_eta)

  if (myrank == 0) then
    print*,'initial kernels:'
    print*,'  bulk min/max : ',min_bulk,max_bulk
    print*,'  betav min/max: ',min_vsv,max_vsv
    print*,'  betah min/max: ',min_vsh,max_vsh
    print*,'  eta min/max  : ',min_eta,max_eta
    print*
  endif

end subroutine read_kernels_tiso

