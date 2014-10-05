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

! sum_kernels_globe
!
! this program can be used for event kernel summation,
! where it sums up transverse isotropic kernel files:
!
!   - proc***_reg1_bulk_c_kernel.bin
!   - proc***_reg1_bulk_betav_kernel.bin
!   - proc***_reg1_bulk_betah_kernel.bin
!   - proc***_reg1_eta_kernel.bin
!
! input file: kernels_run.globe
!   lists all event kernel directories which should be summed together
!
! input directory:  INPUT_KERNELS/
!    contains links to all event kernel directories (listed in "kernel_list.globe")
!
! output directory: OUTPUT_SUM/
!    the resulting kernel files will be stored in this directory

program sum_kernels_globe

  use constants,only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NX_BATHY,NY_BATHY,IIN

  implicit none

  include 'OUTPUT_FILES/values_from_mesher.h'

! ======================================================
  ! USER PARAMETERS

  ! by default, this algorithm uses transverse isotropic (bulk,bulk_betav,bulk_betah,eta) kernels to sum up
  ! if you prefer using isotropic kernels, set flags below accordingly

  ! if you prefer using isotropic kernels (bulk,bulk_beta,rho) kernels, set this flag to true
  logical, parameter :: USE_ISO_KERNELS = .false.

  ! if you prefer isotropic  (alpha,beta,rho) kernels, set this flag to true
  logical, parameter :: USE_ALPHA_BETA_RHO = .false.


! ======================================================

  character(len=150) :: kernel_file_list, kernel_list(1000), sline, kernel_name
  integer :: nker, myrank, sizeprocs
  integer :: ios

  ! ============ program starts here =====================
  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  if(myrank==0) write(*,*) 'reading kernel list: '

  ! default values
  kernel_file_list='kernels_run.globe'

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
     kernel_list(nker) = sline
  enddo
  close(IIN)
  if( myrank == 0 ) then
    write(*,*) '  ',nker,' events'
    write(*,*)
  endif

  ! sums up kernels
  if( USE_ISO_KERNELS ) then

    !  isotropic kernels
    kernel_name = 'reg1_bulk_c_kernel'
    call sum_kernel(kernel_name,kernel_list,nker,myrank)

    kernel_name = 'reg1_bulk_beta_kernel'
    call sum_kernel(kernel_name,kernel_list,nker,myrank)

    kernel_name = 'reg1_rho_kernel'
    call sum_kernel(kernel_name,kernel_list,nker,myrank)

  else if( USE_ALPHA_BETA_RHO ) then

    ! isotropic kernels

    kernel_name = 'reg1_alpha_kernel'
    call sum_kernel(kernel_name,kernel_list,nker,myrank)

    kernel_name = 'reg1_beta_kernel'
    call sum_kernel(kernel_name,kernel_list,nker,myrank)

    kernel_name = 'reg1_rho_kernel'
    call sum_kernel(kernel_name,kernel_list,nker,myrank)

  else

    ! transverse isotropic kernels

    kernel_name = 'reg1_bulk_c_kernel'
    call sum_kernel(kernel_name,kernel_list,nker,myrank)

    kernel_name = 'reg1_bulk_betav_kernel'
    call sum_kernel(kernel_name,kernel_list,nker,myrank)

    kernel_name = 'reg1_bulk_betah_kernel'
    call sum_kernel(kernel_name,kernel_list,nker,myrank)

    kernel_name = 'reg1_eta_kernel'
    call sum_kernel(kernel_name,kernel_list,nker,myrank)

  endif

  if(myrank==0) write(*,*) 'done writing all kernels'

  ! stop all the MPI processes, and exit
  call finalize_mpi()

end program sum_kernels_globe

!
!-------------------------------------------------------------------------------------------------
!

subroutine sum_kernel(kernel_name,kernel_list,nker,myrank)

  use constants,only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NX_BATHY,NY_BATHY,IIN,IOUT

  implicit none

  include 'OUTPUT_FILES/values_from_mesher.h'

  character(len=150) :: kernel_name,kernel_list(1000)
  integer :: nker,myrank

  ! local parameters
  character(len=150) :: k_file
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: &
    kernel_crust_mantle,total_kernel
  double precision :: norm,norm_sum
  integer :: iker,ios

  total_kernel = 0._CUSTOM_REAL
  do iker = 1, nker
    ! user output
    if(myrank==0) then
    write(*,*) 'reading in event kernel for: ',trim(kernel_name)
    write(*,*) '    ',iker, ' out of ', nker
    endif

    ! sensitivity kernel
    kernel_crust_mantle = 0._CUSTOM_REAL
    write(k_file,'(a,i6.6,a)') 'INPUT_KERNELS/'//trim(kernel_list(iker)) &
                          //'/proc',myrank,'_'//trim(kernel_name)//'.bin'

    open(IIN,file=trim(k_file),status='old',form='unformatted',action='read',iostat=ios)
    if( ios /= 0 ) then
     write(*,*) '  kernel not found:',trim(k_file)
     cycle
    endif
    read(IIN) kernel_crust_mantle
    close(IIN)

    ! outputs norm of kernel
    norm = sum( kernel_crust_mantle * kernel_crust_mantle )
    call sum_all_dp(norm,norm_sum)
    if( myrank == 0 ) then
    print*,'  norm kernel: ',sqrt(norm_sum)
    print*
    endif

    ! sums all kernels from each event
    total_kernel(:,:,:,1:NSPEC_CRUST_MANTLE) = &
        total_kernel(:,:,:,1:NSPEC_CRUST_MANTLE) &
        + kernel_crust_mantle(:,:,:,1:NSPEC_CRUST_MANTLE)
  enddo

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

end subroutine sum_kernel

