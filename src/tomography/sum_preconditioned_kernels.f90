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

! sum_preconditioned_kernels_globe
!
! this program can be used for event kernel summation,
! where it sums up transverse isotropic kernel files:
!
!   - proc***_reg1_bulk_c_kernel.bin
!   - proc***_reg1_bulk_betav_kernel.bin
!   - proc***_reg1_bulk_betah_kernel.bin
!   - proc***_reg1_eta_kernel.bin
!
! this version uses the approximate Hessian stored as
!   - proc***_reg1_hess_kernel.bin
!          to precondition the transverse isotropic kernel files
!
! input file: kernels_run.globe
!   lists all event kernel directories which should be summed together
!
! input directory:  INPUT_KERNELS/
!    contains links to all event kernel directories (listed in "kernel_list.globe")
!
! output directory: OUTPUT_SUM/
!    the resulting kernel files will be stored in this directory


program sum_preconditioned_kernels_globe

  use tomography_par

  implicit none

  character(len=150) :: kernel_list(MAX_NUM_NODES), sline, kernel_name
  integer :: nker, sizeprocs
  integer :: ios

  ! ============ program starts here =====================
  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  if(myrank==0) then
    write(*,*) 'sum_preconditioned_kernels_globe:'
    write(*,*)
    write(*,*) 'reading kernel list: '
  endif

  ! reads in event list
  nker=0
  open(unit = IIN, file = trim(kernel_file_list), status = 'old',iostat = ios)
  if (ios /= 0) then
     print *,'Error opening ',trim(kernel_file_list),myrank
     stop 1
  endif
  do while (1 == 1)
     read(IIN,'(a)',iostat=ios) sline
     if (ios /= 0) exit
     nker = nker+1
     if (nker > MAX_NUM_NODES) stop 'Error number of kernels exceeds MAX_NUM_NODES'
     kernel_list(nker) = sline
  enddo
  close(IIN)
  if( myrank == 0 ) then
    write(*,*) '  ',nker,' events'
    write(*,*)
  endif

  ! synchronizes
  call synchronize_all()

  ! sums up kernels
  if( USE_ISO_KERNELS ) then

    !  isotropic kernels
    if( myrank == 0 ) write(*,*) 'isotropic kernels: bulk_c, bulk_beta, rho'

    kernel_name = 'reg1_bulk_c_kernel'
    call sum_kernel_pre(kernel_name,kernel_list,nker)

    kernel_name = 'reg1_bulk_beta_kernel'
    call sum_kernel_pre(kernel_name,kernel_list,nker)

    kernel_name = 'reg1_rho_kernel'
    call sum_kernel_pre(kernel_name,kernel_list,nker)

  else if( USE_ALPHA_BETA_RHO ) then

    ! isotropic kernels
    if( myrank == 0 ) write(*,*) 'isotropic kernels: alpha, beta, rho'

    kernel_name = 'reg1_alpha_kernel'
    call sum_kernel_pre(kernel_name,kernel_list,nker)

    kernel_name = 'reg1_beta_kernel'
    call sum_kernel_pre(kernel_name,kernel_list,nker)

    kernel_name = 'reg1_rho_kernel'
    call sum_kernel_pre(kernel_name,kernel_list,nker)

  else

    ! transverse isotropic kernels
    if( myrank == 0 ) write(*,*) 'transverse isotropic kernels: bulk_c, bulk_betav, bulk_betah,eta'

    kernel_name = 'reg1_bulk_c_kernel'
    call sum_kernel_pre(kernel_name,kernel_list,nker)

    kernel_name = 'reg1_bulk_betav_kernel'
    call sum_kernel_pre(kernel_name,kernel_list,nker)

    kernel_name = 'reg1_bulk_betah_kernel'
    call sum_kernel_pre(kernel_name,kernel_list,nker)

    kernel_name = 'reg1_eta_kernel'
    call sum_kernel_pre(kernel_name,kernel_list,nker)

  endif

  if(myrank==0) write(*,*) 'done writing all kernels'

  ! stop all the MPI processes, and exit
  call finalize_mpi()

end program sum_preconditioned_kernels_globe

!
!-------------------------------------------------------------------------------------------------
!

subroutine sum_kernel_pre(kernel_name,kernel_list,nker)

  use tomography_par

  implicit none

  character(len=150) :: kernel_name,kernel_list(MAX_NUM_NODES)
  integer :: nker

  ! local parameters
  character(len=150) :: k_file
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: kernel,hess,total_kernel
  double precision :: norm,norm_sum
  integer :: iker,ios
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: total_hess,mask_source

  ! initializes arrays
  if( USE_HESS_SUM ) then
    allocate( total_hess(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) )
    total_hess(:,:,:,:) = 0.0_CUSTOM_REAL
  endif

  if( USE_SOURCE_MASK ) then
    allocate( mask_source(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) )
    mask_source(:,:,:,:) = 1.0_CUSTOM_REAL
  endif

  ! loops over all event kernels
  total_kernel = 0._CUSTOM_REAL
  do iker = 1, nker
    ! user output
    if(myrank==0) then
      write(*,*) 'reading in event kernel for: ',trim(kernel_name)
      write(*,*) '  and preconditioner: ','reg1_hess_kernel'
      write(*,*) '    ',iker, ' out of ', nker
    endif

    ! sensitivity kernel / frechet derivative
    kernel = 0._CUSTOM_REAL
    write(k_file,'(a,i6.6,a)') 'INPUT_KERNELS/'//trim(kernel_list(iker)) &
                          //'/proc',myrank,'_'//trim(kernel_name)//'.bin'

    open(IIN,file=trim(k_file),status='old',form='unformatted',action='read',iostat=ios)
    if( ios /= 0 ) then
      write(*,*) '  kernel not found:',trim(k_file)
      cycle
    endif
    read(IIN) kernel
    close(IIN)

    ! outputs norm of kernel
    norm = sum( kernel * kernel )
    call sum_all_dp(norm,norm_sum)
    if( myrank == 0 ) then
      print*,'  norm kernel: ',sqrt(norm_sum)
    endif

    ! approximate Hessian
    hess = 0._CUSTOM_REAL
    write(k_file,'(a,i6.6,a)') 'INPUT_KERNELS/'//trim(kernel_list(iker)) &
                          //'/proc',myrank,'_reg1_hess_kernel.bin'

    open(IIN,file=trim(k_file),status='old',form='unformatted',action='read',iostat=ios)
    if( ios /= 0 ) then
      write(*,*) '  hessian kernel not found: ',trim(k_file)
      stop 'Error hess_kernel.bin files not found'
    endif

    read(IIN) hess
    close(IIN)

    ! outputs norm of preconditioner
    norm = sum( hess * hess )
    call sum_all_dp(norm,norm_sum)
    if( myrank == 0 ) then
      print*,'  norm preconditioner: ',sqrt(norm_sum)
      print*
    endif

    ! note: we take absolute values for hessian (as proposed by Yang)
    hess = abs(hess)

    ! source mask
    if( USE_SOURCE_MASK ) then
      ! reads in mask
      write(k_file,'(a,i6.6,a)') 'INPUT_KERNELS/'//trim(kernel_list(iker)) &
                            //'/proc',myrank,'_reg1_mask_source.bin'
      open(IIN,file=trim(k_file),status='old',form='unformatted',action='read',iostat=ios)
      if( ios /= 0 ) then
        write(*,*) '  file not found:',trim(k_file)
        cycle
      endif
      read(IIN) mask_source
      close(IIN)

      ! masks source elements
      kernel = kernel * mask_source
    endif

    ! precondition
    if( USE_HESS_SUM ) then

      ! sums up hessians first
      total_hess = total_hess + hess

    else

      ! inverts hessian
      call invert_hess( hess )

      ! preconditions each event kernel with its hessian
      kernel = kernel * hess

    endif

    ! sums all kernels from each event
    total_kernel = total_kernel + kernel

  enddo

  ! preconditions summed kernels with summed hessians
  if( USE_HESS_SUM ) then

      ! inverts hessian matrix
      call invert_hess( total_hess )

      ! preconditions kernel
      total_kernel = total_kernel * total_hess

  endif

  ! stores summed kernels
  if(myrank==0) write(*,*) 'writing out summed kernel for: ',trim(kernel_name)

  write(k_file,'(a,i6.6,a)') 'OUTPUT_SUM/proc',myrank,'_'//trim(kernel_name)//'.bin'

  open(IOUT,file=trim(k_file),form='unformatted',status='unknown',action='write',iostat=ios)
  if( ios /= 0 ) then
    write(*,*) 'Error kernel not written:',trim(k_file)
    stop 'Error kernel write'
  endif
  write(IOUT) total_kernel
  close(IOUT)

  if(myrank==0) write(*,*)

  ! frees memory
  if( USE_HESS_SUM ) deallocate(total_hess)
  if( USE_SOURCE_MASK ) deallocate(mask_source)

end subroutine sum_kernel_pre

!
!-------------------------------------------------------------------------------------------------
!

subroutine invert_hess( hess_matrix )

! inverts the hessian matrix
! the approximate hessian is only defined for diagonal elements: like
! H_nn = \frac{ \partial^2 \chi }{ \partial \rho_n \partial \rho_n }
! on all GLL points, which are indexed (i,j,k,ispec)

  use tomography_par

  implicit none

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: hess_matrix

  ! local parameters
  real(kind=CUSTOM_REAL) :: maxh,maxh_all

  ! maximum value of hessian
  maxh = maxval( abs(hess_matrix) )

  ! determines maximum from all slices on master
  call max_allreduce_cr(maxh,maxh_all)

  ! user output
  if( myrank == 0 ) then
    print*
    print*,'hessian maximum: ',maxh_all
    print*
  endif

  ! normalizes hessian
  if( maxh_all < 1.e-18 ) then
    ! hessian is zero, re-initializes
    hess_matrix = 1.0_CUSTOM_REAL
    !call exit_mpi(myrank,'Error hessian too small')
  else
    ! since hessian has absolute values, this scales between [0,1]
    hess_matrix = hess_matrix / maxh_all
  endif


  ! inverts hessian values
  where( abs(hess_matrix(:,:,:,:)) > THRESHOLD_HESS )
    hess_matrix = 1.0_CUSTOM_REAL / hess_matrix
  elsewhere
    hess_matrix = 1.0_CUSTOM_REAL / THRESHOLD_HESS
  endwhere

  ! rescales hessian
  !hess_matrix = hess_matrix * maxh_all

end subroutine invert_hess
