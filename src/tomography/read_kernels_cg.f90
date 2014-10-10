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


subroutine read_kernels_cg_tiso_old()

! reads in smoothed kernels from former iteration in OUTPUT_SUM.old/ : bulk, betav, betah, eta

  use tomography_kernels_tiso_cg

  implicit none
  real(kind=CUSTOM_REAL) :: min_vsv,min_vsh,max_vsv,max_vsh,min_eta,max_eta,min_bulk,max_bulk
  logical:: exist,exist_all,use_old_gradient_all
  integer :: ier
  character(len=150) :: m_file, fname

  ! checks if files are available:
  write(m_file,'(a,i6.6,a)') trim(kernel_old_dir)//'/proc',myrank,'_reg1_bulk_c_kernel_smooth.bin'
  inquire(file=trim(m_file),EXIST=exist)
  if (.not. exist) then
    print*,'Error file does not exist: ',trim(m_file)
    call exit_mpi(myrank,'file not exist')
  endif
  ! makes sure all processes have same flag
  call lor_allreduce_l(exist,exist_all)
  if (.not. exist_all) then
    print*,'old kernels do not exist: ',trim(m_file)
    call exit_mpi(myrank,'flags old model not consistent')
  endif

  ! bulk kernel
  write(m_file,'(a,i6.6,a)') trim(kernel_old_dir)//'/proc',myrank,'_reg1_bulk_c_kernel_smooth.bin'
  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print*,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_bulk_old(:,:,:,1:nspec)
  close(IIN)

  ! betav kernel
  write(m_file,'(a,i6.6,a)') trim(kernel_old_dir)//'/proc',myrank,'_reg1_bulk_betav_kernel_smooth.bin'
  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print*,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_betav_old(:,:,:,1:nspec)
  close(IIN)

  ! betah kernel
  write(m_file,'(a,i6.6,a)') trim(kernel_old_dir)//'/proc',myrank,'_reg1_bulk_betah_kernel_smooth.bin'
  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print*,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_betah_old(:,:,:,1:nspec)
  close(IIN)

  ! eta kernel
  write(m_file,'(a,i6.6,a)') trim(kernel_old_dir)//'/proc',myrank,'_reg1_eta_kernel_smooth.bin'
  open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print*,'Error opening: ',trim(m_file)
    call exit_mpi(myrank,'file not found')
  endif
  read(IIN) kernel_eta_old(:,:,:,1:nspec)
  close(IIN)

  ! statistics
  call min_all_cr(minval(kernel_bulk_old),min_bulk)
  call max_all_cr(maxval(kernel_bulk_old),max_bulk)

  call min_all_cr(minval(kernel_betah_old),min_vsh)
  call max_all_cr(maxval(kernel_betah_old),max_vsh)

  call min_all_cr(minval(kernel_betav_old),min_vsv)
  call max_all_cr(maxval(kernel_betav_old),max_vsv)

  call min_all_cr(minval(kernel_eta_old),min_eta)
  call max_all_cr(maxval(kernel_eta_old),max_eta)

  if (myrank == 0) then
    print*,'old kernels:'
    print*,'  bulk min/max : ',min_bulk,max_bulk
    print*,'  betav min/max: ',min_vsv,max_vsv
    print*,'  betah min/max: ',min_vsh,max_vsh
    print*,'  eta min/max  : ',min_eta,max_eta
    print*
  endif

  ! reads in old gradient directions (phi_(n-1))
  USE_OLD_GRADIENT = .true.

  ! checks if files are available:
  write(m_file,'(a,i6.6,a)') trim(kernel_old_dir)//'/proc',myrank,'_reg1_dbulk_c.bin'
  inquire(file=trim(m_file),EXIST=exist)
  if (.not. exist) then
    print*,'old kernel updates do not exist: ',trim(m_file)
    USE_OLD_GRADIENT = .false.
  endif
  ! makes sure all processes have same flag
  call lor_allreduce_l(exist,exist_all)
  if (.not. exist_all) then
    if (myrank == 0) print*,'old kernel updates do not exist for all: ',trim(m_file)
    call exit_mpi(myrank,'flags old model not consistent')
  endif

  ! makes sure all processes have same flag
  use_old_gradient_all = .false.
  call synchronize_all()

  call lor_allreduce_l(USE_OLD_GRADIENT,use_old_gradient_all)
  if (.not. use_old_gradient_all) then
    print*,'old kernel updates exists, not consistent for all: ',trim(m_file)
    call exit_mpi(myrank,'flags old model not consistent')
  endif

  ! reads in old gradients
  if (USE_OLD_GRADIENT) then
    ! bulk kernel
    fname = 'dbulk_c'
    write(m_file,'(a,i6.6,a)') trim(kernel_old_dir)//'/proc',myrank,'_reg1_'//trim(fname)//'.bin'
    open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print*,'Error opening: ',trim(m_file)
      call exit_mpi(myrank,'file not found')
    endif
    read(IIN) model_dbulk_old(:,:,:,1:nspec)
    close(IIN)

    ! betav kernel
    fname = 'dbetav'
    write(m_file,'(a,i6.6,a)') trim(kernel_old_dir)//'/proc',myrank,'_reg1_'//trim(fname)//'.bin'
    open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print*,'Error opening: ',trim(m_file)
      call exit_mpi(myrank,'file not found')
    endif
    read(IIN) model_dbetav_old(:,:,:,1:nspec)
    close(IIN)

    ! betah kernel
    fname = 'dbetah'
    write(m_file,'(a,i6.6,a)') trim(kernel_old_dir)//'/proc',myrank,'_reg1_'//trim(fname)//'.bin'
    open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print*,'Error opening: ',trim(m_file)
      call exit_mpi(myrank,'file not found')
    endif
    read(IIN) model_dbetah_old(:,:,:,1:nspec)
    close(IIN)

    ! eta kernel
    fname = 'deta'
    write(m_file,'(a,i6.6,a)') trim(kernel_old_dir)//'/proc',myrank,'_reg1_'//trim(fname)//'.bin'
    open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print*,'Error opening: ',trim(m_file)
      call exit_mpi(myrank,'file not found')
    endif
    read(IIN) model_deta_old(:,:,:,1:nspec)
    close(IIN)


    ! statistics
    call min_all_cr(minval(model_dbulk_old),min_bulk)
    call max_all_cr(maxval(model_dbulk_old),max_bulk)

    call min_all_cr(minval(model_dbetah_old),min_vsh)
    call max_all_cr(maxval(model_dbetah_old),max_vsh)

    call min_all_cr(minval(model_dbetav_old),min_vsv)
    call max_all_cr(maxval(model_dbetav_old),max_vsv)

    call min_all_cr(minval(model_deta_old),min_eta)
    call max_all_cr(maxval(model_deta_old),max_eta)

    if (myrank == 0) then
      print*,'old kernel updates:'
      print*,'  bulk min/max : ',min_bulk,max_bulk
      print*,'  betav min/max: ',min_vsv,max_vsv
      print*,'  betah min/max: ',min_vsh,max_vsh
      print*,'  eta min/max  : ',min_eta,max_eta
      print*
    endif
  endif ! USE_OLD_GRADIENT

end subroutine read_kernels_cg_tiso_old

